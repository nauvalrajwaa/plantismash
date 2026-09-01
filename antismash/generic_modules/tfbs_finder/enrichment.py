"""
CRE (Cis-Regulatory Element) and TFBS cluster-level motif enrichment analysis.

Computes TF family-level enrichment across candidate BGC promoter windows (scoped
to quorum genes from the attribution layer, falling back to all cluster CDS).
Background rates are derived once per record from non-cluster promoter windows.
"""

from __future__ import annotations
from dataclasses import dataclass
import os
import math
import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio.SeqRecord import SeqRecord
from antismash import utils
from .tfbs_detection import (
    _HAS_MOODS,
    PWM_PATH,
    Matrix,
    TFBSHit,
    _merge_intervals,
    _load_matrices_cached,
    _scan_segment_with_pwms,
    build_upstream_caps,
    get_cds_capped_promoter_window,
)

# Path to Ath TF motif information metadata
INFO_PATH = utils.get_full_path(__file__, os.path.join("data", "Ath_TF_binding_motifs_information.txt"))


def load_tf_family_mapping(info_path: str = INFO_PATH) -> Dict[str, str]:
    """
    Load mapping from motif / gene id to TF Family from Ath_TF_binding_motifs_information.txt.
    Key in PWM DB matches Gene_id (e.g. AT1G01060).
    """
    mapping: Dict[str, str] = {}
    if not os.path.exists(info_path):
        logging.warning("TFBS enrichment: info file missing at %s", info_path)
        return mapping

    try:
        with open(info_path, "r", encoding="utf-8") as fh:
            header = fh.readline()
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    gene_id = parts[0].strip()
                    family = parts[1].strip()
                    if gene_id and family:
                        mapping[gene_id] = family
                    if len(parts) >= 3:
                        matrix_id = parts[2].strip()
                        if matrix_id and family:
                            mapping[matrix_id] = family
    except Exception as e:
        logging.warning("TFBS enrichment: failed to load family mapping: %s", e)

    return mapping


def get_cluster_quorum_gene_ids(cluster_feature, record: SeqRecord) -> Set[str]:
    """
    Retrieve quorum / biosynthetic gene IDs for a cluster from the attribution layer or
    sec_met qualifiers. Fallback to all CDS IDs overlapping the cluster.
    """
    quorum_ids: Set[str] = set()

    # 1. Try sec_met qualifiers in the cluster
    for cds in utils.get_cluster_cds_features(cluster_feature, record):
        gid = utils.get_gene_id(cds)
        quals = getattr(cds, "qualifiers", {})
        if "sec_met" in quals:
            for sm in quals["sec_met"]:
                if sm.startswith("Type:") or sm.startswith("NRPS/PKS") or "bitscore" in sm:
                    quorum_ids.add(gid)
                    break

    # 2. If nothing found in sec_met, fallback to all CDS in cluster
    if not quorum_ids:
        for cds in utils.get_cluster_cds_features(cluster_feature, record):
            quorum_ids.add(utils.get_gene_id(cds))

    return quorum_ids


def log_poisson_pmf(k: int, lambd: float) -> float:
    """Compute ln(P(X = k | lambd)) using math.lgamma."""
    if lambd <= 0:
        return 0.0 if k == 0 else -float("inf")
    if k < 0:
        return -float("inf")
    return k * math.log(lambd) - lambd - math.lgamma(k + 1)


def poisson_upper_tail(k: int, lambd: float) -> float:
    """
    Compute P(X >= k | lambd) in pure Python.

    For k <= lambd the CDF subtraction is numerically safe (the lower CDF is
    bounded well away from 1). For k > lambd a direct upper-tail log-sum-exp is
    used instead of 1 - CDF, which would suffer catastrophic cancellation and
    collapse extreme tails to exactly 0.0.
    """
    if k <= 0:
        return 1.0
    if lambd <= 0:
        return 0.0

    if k <= lambd:
        cdf_lower = 0.0
        for i in range(k):
            cdf_lower += math.exp(log_poisson_pmf(i, lambd))
        return max(0.0, min(1.0, 1.0 - cdf_lower))

    # Direct upper-tail summation with truncation
    log_terms = []
    curr = k
    while True:
        lp = log_poisson_pmf(curr, lambd)
        log_terms.append(lp)
        # stop if contribution is negligible compared to peak or curr far past lambd
        if curr > lambd and (lp < log_terms[0] - 30.0 or curr > k + 2000):
            break
        curr += 1

    max_lp = max(log_terms)
    p = sum(math.exp(lp - max_lp) for lp in log_terms) * math.exp(max_lp)
    return max(0.0, min(1.0, p))


def benjamini_hochberg_fdr(p_values: List[float]) -> List[float]:
    """
    Benjamini-Hochberg FDR correction over a list of p-values.
    Returns q-values in original order.
    """
    n = len(p_values)
    if n == 0:
        return []

    sorted_indices = sorted(range(n), key=lambda i: p_values[i])
    q_values = [1.0] * n

    min_q = 1.0
    for rank_minus_1 in range(n - 1, -1, -1):
        idx = sorted_indices[rank_minus_1]
        rank = rank_minus_1 + 1
        q = (p_values[idx] * n) / rank
        if q < min_q:
            min_q = q
        q_values[idx] = min(1.0, max(0.0, min_q))

    return q_values


def compute_record_background_rates(record: SeqRecord,
                                    upstream_bp: int,
                                    pvalue: float,
                                    matrices: List[Matrix]) -> Dict[str, float]:
    """
    Precompute per-motif background hit rates (hits/kb) across all non-cluster promoter windows.
    Runs once per record.
    """
    if not _HAS_MOODS or not matrices:
        return {}

    seqlen = len(record.seq)
    clusters = utils.get_sorted_cluster_features(record)
    all_cds = list(utils.get_cds_features(record))

    # Exclude CDS overlapping any cluster
    non_bgc_cds = []
    for cds in all_cds:
        c_start = int(cds.location.start)
        c_end = int(cds.location.end)
        in_cluster = False
        for cl in clusters:
            if c_end > int(cl.location.start) and c_start < int(cl.location.end):
                in_cluster = True
                break
        if not in_cluster:
            non_bgc_cds.append(cds)

    if not non_bgc_cds:
        return {}

    caps = build_upstream_caps(all_cds)

    raw_windows = []
    for cds in non_bgc_cds:
        w = get_cds_capped_promoter_window(cds, caps, upstream_bp, seqlen)
        if w is not None:
            raw_windows.append(w)

    merged_windows = _merge_intervals(raw_windows)
    total_kb = sum(b - a + 1 for a, b in merged_windows) / 1000.0
    if total_kb <= 0:
        return {}

    motif_hits: Dict[str, int] = {m.name: 0 for m in matrices}
    for a, b in merged_windows:
        seg_hits = _scan_segment_with_pwms(record, a, b, matrices, pvalue)
        for mat_idx, _, _, _ in seg_hits:
            mname = matrices[mat_idx].name
            motif_hits[mname] = motif_hits.get(mname, 0) + 1

    rates = {mname: count / total_kb for mname, count in motif_hits.items()}
    return rates


@dataclass
class CREEnrichmentResult:
    record: str
    cluster: str
    product: str
    n_quorum_genes: int
    promoter_kb_scanned: float
    tf_family: str
    n_motifs: int
    hits_obs: int
    hits_exp: float
    fold_enrichment: float
    p_value: float
    q_value: float
    coherence_fraction: float
    top_motifs: str
    max_confidence: str

    def to_tsv_row(self) -> List[str]:
        fold_str = "inf" if (self.hits_exp == 0 and self.hits_obs > 0) else f"{self.fold_enrichment:.2f}"
        return [
            self.record,
            str(self.cluster),
            self.product,
            str(self.n_quorum_genes),
            f"{self.promoter_kb_scanned:.3f}",
            self.tf_family,
            str(self.n_motifs),
            str(self.hits_obs),
            f"{self.hits_exp:.2f}",
            fold_str,
            f"{self.p_value:.4e}",
            f"{self.q_value:.4e}",
            f"{self.coherence_fraction:.2f}",
            self.top_motifs,
            self.max_confidence,
        ]


def run_cre_enrichment(record: SeqRecord,
                       options,
                       hits: List[TFBSHit],
                       upstream_bp: int,
                       pvalue: float,
                       matrix_path: str = PWM_PATH) -> List[CREEnrichmentResult]:
    """
    Run CRE motif enrichment analysis for all clusters in record.
    Returns list of CREEnrichmentResult objects.
    """
    if not _HAS_MOODS:
        return []

    matrices = _load_matrices_cached(matrix_path)
    if not matrices:
        return []

    family_map = load_tf_family_mapping()
    motif_to_family = {m.name: family_map.get(m.name, "Unknown") for m in matrices}
    all_families: Set[str] = set(motif_to_family.values())
    family_motifs: Dict[str, List[str]] = {f: [] for f in all_families}
    for mname, f in motif_to_family.items():
        family_motifs[f].append(mname)

    # 1. Background rates (precomputed once per record)
    bg_rates = compute_record_background_rates(record, upstream_bp, pvalue, matrices)

    seqlen = len(record.seq)
    clusters = utils.get_sorted_cluster_features(record)
    all_cds = list(utils.get_cds_features(record))

    caps = build_upstream_caps(all_cds)

    cluster_tests: List[Dict[str, Any]] = []

    for cluster in clusters:
        cnum = utils.get_cluster_number(cluster)
        cstart = int(cluster.location.start)
        cend = int(cluster.location.end)
        ctype = utils.get_cluster_type(cluster) if hasattr(utils, 'get_cluster_type') else "-"
        if not ctype or ctype == "-":
            ctype = cluster.qualifiers.get("product", ["-"])[0]

        quorum_gids = get_cluster_quorum_gene_ids(cluster, record)
        quorum_cds_list = []
        quorum_windows = []
        cds_window_map = {}

        cluster_cds_feats = [
            f for f in all_cds
            if int(f.location.end) > cstart and int(f.location.start) < cend
        ]

        for cds in cluster_cds_feats:
            gid = utils.get_gene_id(cds)
            if gid in quorum_gids:
                w = get_cds_capped_promoter_window(cds, caps, upstream_bp, seqlen)
                if w is not None:
                    # Clip window to cluster boundaries
                    a_clip = max(w[0], cstart)
                    b_clip = min(w[1], cend - 1)
                    if b_clip >= a_clip:
                        w_clipped = (a_clip, b_clip)
                        quorum_cds_list.append((gid, cds, w_clipped))
                        quorum_windows.append(w_clipped)
                        cds_window_map[gid] = w_clipped

        if not quorum_cds_list:
            continue

        merged_q_windows = _merge_intervals(quorum_windows)
        kb_scanned = sum(b - a + 1 for a, b in merged_q_windows) / 1000.0
        n_quorum = len(quorum_cds_list)

        # Map observed hits to quorum CDS windows
        quorum_hits: List[Tuple[str, TFBSHit]] = []  # (gene_id, hit)
        for gid, _, (wa, wb) in quorum_cds_list:
            for h in hits:
                mlen = len(h.consensus) if h.consensus else 1
                if wa <= h.start <= max(wa, wb - mlen):
                    quorum_hits.append((gid, h))

        # Evaluate per family
        for fam, mnames in family_motifs.items():
            fam_hits = [h for gid, h in quorum_hits if h.name in mnames]
            obs = len(fam_hits)

            # Expected hits = sum_{m in fam} rate_m * kb_scanned
            fam_bg_rate = sum(bg_rates.get(m, 0.0) for m in mnames)
            exp = fam_bg_rate * kb_scanned

            # Skip families with zero genome-wide hits AND zero cluster hits
            if obs == 0 and exp == 0:
                continue

            # Coherence
            genes_with_hit = set(gid for gid, h in quorum_hits if h.name in mnames)
            coherence = len(genes_with_hit) / float(n_quorum) if n_quorum > 0 else 0.0

            # Top motifs by score
            fam_hits_sorted = sorted(fam_hits, key=lambda h: h.score, reverse=True)
            top_m_names = []
            seen_m = set()
            for h in fam_hits_sorted:
                if h.name not in seen_m:
                    seen_m.add(h.name)
                    top_m_names.append(h.name)
                if len(top_m_names) >= 3:
                    break
            top_motifs_str = ";".join(top_m_names) if top_m_names else "-"

            # Max confidence
            max_conf = "none"
            if fam_hits:
                conf_order = {"strong": 3, "medium": 2, "weak": 1}
                best_h = max(fam_hits, key=lambda h: conf_order.get(str(h.confidence).lower(), 0))
                max_conf = str(best_h.confidence).capitalize()

            # Fold enrichment
            fold = obs / exp if exp > 0 else (float("inf") if obs > 0 else 1.0)

            # Poisson upper tail
            p_val = poisson_upper_tail(obs, exp)

            cluster_tests.append({
                "record": record.id,
                "cluster": str(cnum),
                "product": ctype,
                "n_quorum_genes": n_quorum,
                "promoter_kb_scanned": kb_scanned,
                "tf_family": fam,
                "n_motifs": len(mnames),
                "hits_obs": obs,
                "hits_exp": exp,
                "fold_enrichment": fold,
                "p_value": p_val,
                "coherence_fraction": coherence,
                "top_motifs": top_motifs_str,
                "max_confidence": max_conf,
            })

    if not cluster_tests:
        return []

    # 2. Multiple testing correction (BH-FDR across all cluster x family tests)
    all_p_values = [t["p_value"] for t in cluster_tests]
    all_q_values = benjamini_hochberg_fdr(all_p_values)

    results: List[CREEnrichmentResult] = []
    for i, t in enumerate(cluster_tests):
        results.append(
            CREEnrichmentResult(
                record=t["record"],
                cluster=t["cluster"],
                product=t["product"],
                n_quorum_genes=t["n_quorum_genes"],
                promoter_kb_scanned=t["promoter_kb_scanned"],
                tf_family=t["tf_family"],
                n_motifs=t["n_motifs"],
                hits_obs=t["hits_obs"],
                hits_exp=t["hits_exp"],
                fold_enrichment=t["fold_enrichment"],
                p_value=t["p_value"],
                q_value=all_q_values[i],
                coherence_fraction=t["coherence_fraction"],
                top_motifs=t["top_motifs"],
                max_confidence=t["max_confidence"],
            )
        )

    # Sort by cluster index then q_value ascending
    def sort_key(r: CREEnrichmentResult):
        try:
            c_int = int(r.cluster)
        except ValueError:
            c_int = 999999
        return (c_int, r.q_value)

    results.sort(key=sort_key)
    return results


CRE_TSV_HEADER = [
    "record", "cluster", "product", "n_quorum_genes", "promoter_kb_scanned",
    "tf_family", "n_motifs", "hits_obs", "hits_exp", "fold_enrichment",
    "p_value", "q_value", "coherence_fraction", "top_motifs", "max_confidence",
]


def write_cre_enrichment_tsv(results: List[CREEnrichmentResult], outdir: str) -> str:
    """
    Write CRE enrichment results to <output>/tfbs/cre_enrichment.tsv.
    Always writes header even when results list is empty.
    """
    tfbs_dir = os.path.join(outdir, "tfbs")
    os.makedirs(tfbs_dir, exist_ok=True)
    tsv_path = os.path.join(tfbs_dir, "cre_enrichment.tsv")

    try:
        with open(tsv_path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(CRE_TSV_HEADER) + "\n")
            for r in results:
                fh.write("\t".join(r.to_tsv_row()) + "\n")
        logging.info("TFBS: written CRE enrichment TSV to %s (%d rows)", tsv_path, len(results))
    except Exception as e:
        logging.warning("TFBS: failed to write CRE enrichment TSV: %s", e)

    return tsv_path
