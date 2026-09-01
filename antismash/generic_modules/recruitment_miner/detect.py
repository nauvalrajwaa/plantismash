"""
Pure logic and detection functions for recruitment_miner (Feature B).

Identifies primary-metabolism / essential-gene paralog recruitments inside BGCs:
  - CDS vs Essential DB alignment filtering
  - Quorum-exclusion check (genes outside the biosynthetic core)
  - Genome duplication check (must have paralogs outside the candidate BGC)
  - Pseudogenization / corruption tagging (qualifiers, internal stops, truncation)
"""

from __future__ import annotations
from dataclasses import dataclass
import csv
import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from antismash import utils


# Default thresholds for recruitment mining
DEFAULT_EVALUE = 1e-10
DEFAULT_PIDENT = 40.0
DEFAULT_QCOV = 60.0


@dataclass
class EssentialHit:
    record_id: str
    cluster_idx: int
    product: str
    gene_id: str
    agi_hit: str
    pident: float
    qcov: float
    evalue: float
    source: str
    family: str
    dup_outside: str          # "yes" | "no"
    copies_outside: int
    corrupted: str            # "yes" | "no"
    corruption_reason: str
    quorum_excluded: str      # "yes" | "no"


@dataclass
class ClusterRecruitmentFlag:
    record_id: str
    cluster_idx: int
    product: str
    n_candidates: int
    n_corrupted: int
    gene_ids: List[str]
    bonus_signal: str         # e.g. "2 paralogs (1 corrupted)"


@dataclass
class AlignmentMatch:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    evalue: float
    qlen: int
    slen: int

    @property
    def qcov(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return (self.length / float(self.qlen)) * 100.0


def locate_diamond_binary() -> Optional[str]:
    """Find diamond executable on PATH or in conda environment."""
    binary = utils.locate_executable("diamond")
    if binary:
        return binary

    # Check common fallback locations
    candidates = [
        shutil.which("diamond"),
        os.path.expanduser("~/miniforge3/envs/plantismash/bin/diamond"),
        os.path.expanduser("~/miniconda3/envs/plantismash/bin/diamond"),
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


def parse_diamond_tabular(tabular_text: str) -> List[AlignmentMatch]:
    """
    Parse DIAMOND blastp tabular output:
    qseqid sseqid pident length evalue qlen slen
    """
    matches: List[AlignmentMatch] = []
    for line in tabular_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        try:
            qseqid = parts[0].strip()
            sseqid = parts[1].strip()
            pident = float(parts[2])
            length = int(parts[3])
            evalue = float(parts[4])
            qlen = int(parts[5])
            slen = int(parts[6])
            matches.append(
                AlignmentMatch(
                    qseqid=qseqid,
                    sseqid=sseqid,
                    pident=pident,
                    length=length,
                    evalue=evalue,
                    qlen=qlen,
                    slen=slen,
                )
            )
        except (ValueError, IndexError):
            continue
    return matches


def load_reference_metadata(metadata_tsv_path: str) -> Dict[str, Tuple[str, str, str]]:
    """
    Load reference metadata TSV:
    agi_locus -> (source, family, evidence_note)
    """
    meta: Dict[str, Tuple[str, str, str]] = {}
    if not os.path.isfile(metadata_tsv_path):
        return meta

    try:
        with open(metadata_tsv_path, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            for row in reader:
                if not row or len(row) < 3:
                    continue
                locus = row[0].strip().upper()
                source = row[1].strip()
                family = row[2].strip()
                note = row[3].strip() if len(row) > 3 else ""
                meta[locus] = (source, family, note)
    except Exception as e:
        logging.warning("recruitment_miner: failed to read metadata TSV %s: %s", metadata_tsv_path, e)

    return meta


def get_cluster_quorum_gene_ids(cluster_feature: SeqFeature, record: SeqRecord) -> Set[str]:
    """
    Retrieve quorum/biosynthetic gene IDs for a cluster from sec_met qualifiers.
    Matches the proxy used in CRE enrichment.
    """
    quorum_ids: Set[str] = set()
    for cds in utils.get_cluster_cds_features(cluster_feature, record):
        gid = utils.get_gene_id(cds)
        quals = getattr(cds, "qualifiers", {})
        if "sec_met" in quals:
            for sm in quals["sec_met"]:
                if sm.startswith("Type:") or sm.startswith("NRPS/PKS") or "bitscore" in sm:
                    quorum_ids.add(gid)
                    break
    return quorum_ids


def check_corruption(
    feature: SeqFeature,
    db_match: Optional[AlignmentMatch] = None,
    outside_matches: Optional[List[AlignmentMatch]] = None,
) -> Tuple[bool, str]:
    """
    Evaluate pseudogenization / corruption signals on a CDS feature:
      1. /pseudo, /pseudogene, /partial qualifier present
      2. Internal '*' in translation
      3. essential-DB alignment qcov < 90% while paralogs outside have qcov >= 90% (truncated copy)
    Returns (is_corrupted, reason_string).
    """
    quals = getattr(feature, "qualifiers", {})
    reasons = []

    # 1. Qualifiers
    if "pseudo" in quals:
        reasons.append("/pseudo qualifier")
    if "pseudogene" in quals:
        reasons.append("/pseudogene qualifier")
    if "partial" in quals:
        reasons.append("/partial qualifier")

    # 2. Translation internal stop codon
    translations = quals.get("translation", [])
    for trans in translations:
        if isinstance(trans, str):
            # Internal stop: '*' appears before the end of the sequence
            stripped = trans.rstrip("*")
            if "*" in stripped:
                reasons.append("internal stop codon in translation")
                break

    # 3. Truncation relative to genome copies
    if db_match and outside_matches:
        if db_match.qcov < 90.0:
            for om in outside_matches:
                if om.qcov >= 90.0 and om.pident >= 40.0:
                    reasons.append(
                        f"truncated alignment (cluster qcov={db_match.qcov:.1f}% vs outside copy {om.qseqid} qcov={om.qcov:.1f}%)"
                    )
                    break

    if reasons:
        return True, "; ".join(reasons)
    return False, "clean"


def run_diamond_blastp(
    query_faa: str,
    db_path: str,
    threads: int = 1,
    evalue: float = DEFAULT_EVALUE,
) -> List[AlignmentMatch]:
    """
    Run DIAMOND blastp query_faa vs db_path with tabular output:
    qseqid sseqid pident length evalue qlen slen
    """
    diamond_bin = locate_diamond_binary()
    if not diamond_bin:
        logging.warning("recruitment_miner: DIAMOND binary not found; skipping alignment.")
        return []

    # Strip .dmnd extension if supplied for db
    db_base = db_path[:-5] if db_path.endswith(".dmnd") else db_path

    cmd = [
        diamond_bin,
        "blastp",
        "--db",
        db_base,
        "--query",
        query_faa,
        "--threads",
        str(max(1, threads)),
        "--evalue",
        str(evalue),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "evalue",
        "qlen",
        "slen",
        "--max-target-seqs",
        "10000",
    ]

    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return parse_diamond_tabular(res.stdout)
    except subprocess.CalledProcessError as e:
        logging.warning("recruitment_miner: DIAMOND run failed (code %d): %s", e.returncode, e.stderr)
        return []


def analyze_record_recruitment(
    record: SeqRecord,
    options: Any,
    db_dir: str,
    evalue_cutoff: float = DEFAULT_EVALUE,
    pident_cutoff: float = DEFAULT_PIDENT,
    qcov_cutoff: float = DEFAULT_QCOV,
) -> Tuple[List[EssentialHit], List[ClusterRecruitmentFlag]]:
    """
    Perform recruitment mining on a single SeqRecord.
    Returns (list_of_essential_hits, list_of_cluster_flags).
    """
    essential_hits: List[EssentialHit] = []
    cluster_flags: List[ClusterRecruitmentFlag] = []

    # 1. Collect all CDS features with valid translations
    cds_features = utils.get_cds_features(record)
    valid_cds: List[Tuple[str, SeqFeature, str]] = []  # (gene_id, feature, translation)
    for cds in cds_features:
        gid = utils.get_gene_id(cds)
        trans = cds.qualifiers.get("translation", [""])[0]
        if trans and isinstance(trans, str) and len(trans.strip()) > 0:
            valid_cds.append((gid, cds, trans.strip()))

    if not valid_cds:
        return essential_hits, cluster_flags

    # 2. Check essential database presence
    dmnd_path = os.path.join(db_dir, "essential_proteins.dmnd")
    faa_path = os.path.join(db_dir, "essential_proteins.faa")
    meta_path = os.path.join(db_dir, "reference_metadata.tsv")

    if not (os.path.isfile(dmnd_path) or os.path.isfile(faa_path)):
        logging.warning(
            "recruitment_miner: essential database not found at %s. Returning empty results.", db_dir
        )
        return essential_hits, cluster_flags

    # Auto-compile .dmnd if only .faa exists. Inlined diamond makedb call:
    # bash_scripts/ is not packaged into site-packages on server installs
    # (non-editable pip install), so it must never be imported at runtime.
    if not os.path.isfile(dmnd_path) and os.path.isfile(faa_path):
        diamond_bin = locate_diamond_binary()
        if not diamond_bin:
            logging.warning("recruitment_miner: DIAMOND binary not found; cannot build .dmnd from .faa.")
            return essential_hits, cluster_flags
        db_base = dmnd_path[:-5] if dmnd_path.endswith(".dmnd") else dmnd_path
        try:
            subprocess.run(
                [diamond_bin, "makedb", "--in", faa_path, "-d", db_base],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True,
            )
        except subprocess.CalledProcessError as e:
            logging.warning("recruitment_miner: diamond makedb failed (code %d): %s", e.returncode, e.stderr)
            return essential_hits, cluster_flags

    metadata_map = load_reference_metadata(meta_path)
    threads = getattr(options, "cpus", 1)

    # 3. Create temp directory for FASTA & alignment
    temp_dir = tempfile.mkdtemp(prefix="plantismash_recruitment_")
    try:
        record_faa = os.path.join(temp_dir, "record_cds.faa")
        with open(record_faa, "w", encoding="utf-8") as fh:
            for gid, _, trans in valid_cds:
                fh.write(f">{gid}\n{trans}\n")

        # Alignment 1: Record CDS vs Essential DB
        db_matches = run_diamond_blastp(record_faa, dmnd_path, threads=threads, evalue=evalue_cutoff)

        # Build self-DB for duplication check
        self_dmnd = os.path.join(temp_dir, "record_self.dmnd")
        diamond_bin = locate_diamond_binary()
        if diamond_bin:
            cmd_mk = [diamond_bin, "makedb", "--in", record_faa, "-d", self_dmnd]
            subprocess.run(cmd_mk, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Alignment 2: Record CDS vs Self
        self_matches = (
            run_diamond_blastp(record_faa, self_dmnd, threads=threads, evalue=evalue_cutoff)
            if os.path.isfile(self_dmnd + ".dmnd") or os.path.isfile(self_dmnd)
            else []
        )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    # Map best DB hits per query gene: gid -> best AlignmentMatch
    best_db_hits: Dict[str, AlignmentMatch] = {}
    for m in db_matches:
        if m.pident >= pident_cutoff and m.qcov >= qcov_cutoff and m.evalue <= evalue_cutoff:
            if m.qseqid not in best_db_hits or m.evalue < best_db_hits[m.qseqid].evalue:
                best_db_hits[m.qseqid] = m

    # Map self hits per query gene: gid -> list of AlignmentMatch (excluding identity self-hit)
    self_hits_by_gene: Dict[str, List[AlignmentMatch]] = {}
    for sm in self_matches:
        if sm.qseqid != sm.sseqid and sm.pident >= pident_cutoff and sm.qcov >= qcov_cutoff:
            self_hits_by_gene.setdefault(sm.qseqid, []).append(sm)

    # 4. Evaluate clusters in record
    clusters = utils.get_sorted_cluster_features(record)
    feature_by_id = {gid: feat for gid, feat, _ in valid_cds}

    for cluster in clusters:
        cluster_idx = utils.get_cluster_number(cluster)
        product = cluster.qualifiers.get("product", ["unknown"])[0]
        cluster_cds = utils.get_cluster_cds_features(cluster, record)
        cluster_gids = {utils.get_gene_id(c) for c in cluster_cds}
        quorum_gids = get_cluster_quorum_gene_ids(cluster, record)

        cand_in_cluster: List[str] = []
        corrupted_in_cluster: List[str] = []

        for cds in cluster_cds:
            gid = utils.get_gene_id(cds)
            is_quorum = gid in quorum_gids
            quorum_excluded = "no" if is_quorum else "yes"

            db_match = best_db_hits.get(gid)
            if not db_match:
                continue

            # Target AGI info
            target_locus = db_match.sseqid.split()[0].split(".")[0].upper()
            source, family, _ = metadata_map.get(
                target_locus, ("UNKNOWN", "essential_paralog", "")
            )

            # Duplication check outside cluster
            self_hits = self_hits_by_gene.get(gid, [])
            outside_hits = [h for h in self_hits if h.sseqid not in cluster_gids]
            dup_outside = "yes" if len(outside_hits) > 0 else "no"
            copies_outside = len({h.sseqid for h in outside_hits})

            # Corruption check
            is_corrupted, corruption_reason = check_corruption(
                cds, db_match=db_match, outside_matches=outside_hits
            )
            corrupted_str = "yes" if is_corrupted else "no"

            # Check if this qualifies as a recruitment candidate:
            # Must be OUTSIDE quorum AND have duplication outside cluster
            is_candidate = (not is_quorum) and (dup_outside == "yes")
            if is_candidate:
                cand_in_cluster.append(gid)
                if is_corrupted:
                    corrupted_in_cluster.append(gid)

            hit_entry = EssentialHit(
                record_id=record.id,
                cluster_idx=cluster_idx,
                product=product,
                gene_id=gid,
                agi_hit=target_locus,
                pident=db_match.pident,
                qcov=db_match.qcov,
                evalue=db_match.evalue,
                source=source,
                family=family,
                dup_outside=dup_outside,
                copies_outside=copies_outside,
                corrupted=corrupted_str,
                corruption_reason=corruption_reason,
                quorum_excluded=quorum_excluded,
            )
            essential_hits.append(hit_entry)

        # Build cluster flag summary
        n_cand = len(cand_in_cluster)
        n_corr = len(corrupted_in_cluster)
        bonus_sig = f"{n_cand} paralogs ({n_corr} corrupted)" if n_cand > 0 else "0 paralogs (0 corrupted)"

        c_flag = ClusterRecruitmentFlag(
            record_id=record.id,
            cluster_idx=cluster_idx,
            product=product,
            n_candidates=n_cand,
            n_corrupted=n_corr,
            gene_ids=cand_in_cluster,
            bonus_signal=bonus_sig,
        )
        cluster_flags.append(c_flag)

    return essential_hits, cluster_flags
