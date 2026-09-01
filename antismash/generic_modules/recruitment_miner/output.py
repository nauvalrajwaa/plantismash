"""
Output generation and HTML formatting for recruitment_miner (Feature B).

Writes:
  1. <output>/recruitment_miner/essential_hits.tsv
  2. <output>/recruitment_miner/cluster_flags.tsv
Renders cluster HTML detail blocks using PyQuery.
"""

from __future__ import annotations
import csv
import html
import logging
import os
from typing import Any, Dict, List, Optional, Tuple
from pyquery import PyQuery as pq

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from .detect import EssentialHit, ClusterRecruitmentFlag


ESSENTIAL_HITS_HEADER = [
    "record",
    "cluster",
    "product",
    "gene_id",
    "agi_hit",
    "pident",
    "qcov",
    "evalue",
    "source",
    "family",
    "dup_outside",
    "copies_outside",
    "corrupted",
    "corruption_reason",
    "quorum_excluded",
]

CLUSTER_FLAGS_HEADER = [
    "record",
    "cluster",
    "product",
    "n_candidates",
    "n_corrupted",
    "gene_ids",
    "bonus_signal",
]


def _get_output_dir(options: Any) -> str:
    """Return the base results directory."""
    val = getattr(options, "full_outputfolder_path", None)
    if isinstance(val, str) and val.strip():
        return val

    val = getattr(options, "outputfoldername", None)
    if isinstance(val, str) and val.strip():
        return os.path.abspath(val)

    candidates = [
        "outputfolder", "output_dir", "outdir", "output", "output_folder",
        "output_path", "results_dir", "results_path", "result_dir", "work_dir",
    ]
    for attr in candidates:
        v = getattr(options, attr, None)
        if isinstance(v, str) and v.strip():
            return os.path.abspath(v)

    return os.getcwd()


def _ensure_recruitment_dir(options: Any) -> str:
    outdir = _get_output_dir(options)
    target_dir = os.path.join(outdir, "recruitment_miner")
    os.makedirs(target_dir, exist_ok=True)
    return target_dir


def write_recruitment_tsvs(
    hits: List[EssentialHit],
    flags: List[ClusterRecruitmentFlag],
    outdir: str,
) -> Tuple[str, str]:
    """
    Write essential_hits.tsv and cluster_flags.tsv to outdir/recruitment_miner/.
    Always creates files with headers even if empty.
    """
    target_dir = os.path.join(outdir, "recruitment_miner")
    os.makedirs(target_dir, exist_ok=True)

    hits_tsv = os.path.join(target_dir, "essential_hits.tsv")
    flags_tsv = os.path.join(target_dir, "cluster_flags.tsv")

    # Write essential_hits.tsv
    with open(hits_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(ESSENTIAL_HITS_HEADER)
        for h in hits:
            writer.writerow([
                h.record_id,
                h.cluster_idx,
                h.product,
                h.gene_id,
                h.agi_hit,
                f"{h.pident:.1f}",
                f"{h.qcov:.1f}",
                f"{h.evalue:.2e}",
                h.source,
                h.family,
                h.dup_outside,
                h.copies_outside,
                h.corrupted,
                h.corruption_reason,
                h.quorum_excluded,
            ])

    # Write cluster_flags.tsv
    with open(flags_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(CLUSTER_FLAGS_HEADER)
        for f in flags:
            writer.writerow([
                f.record_id,
                f.cluster_idx,
                f.product,
                f.n_candidates,
                f.n_corrupted,
                ";".join(f.gene_ids),
                f.bonus_signal,
            ])

    return hits_tsv, flags_tsv


def generate_details_div(
    cluster_feature: SeqFeature,
    seq_record: SeqRecord,
    options: Any,
    js_domains: Any = None,
    details: Any = None,
) -> Optional[pq]:
    """
    Render HTML section for cluster details view.
    Only renders when this cluster has >=1 candidate paralog.
    """
    try:
        extra_ns = getattr(options, "extrarecord", {}).get(seq_record.id)
        if not extra_ns:
            return details

        extradata = getattr(extra_ns, "extradata", {})
        cluster_flags: List[ClusterRecruitmentFlag] = extradata.get("RecruitmentClusterFlags", [])
        essential_hits: List[EssentialHit] = extradata.get("RecruitmentEssentialHits", [])

        # Find flag for this cluster
        from antismash import utils
        c_num = utils.get_cluster_number(cluster_feature)
        matching_flags = [f for f in cluster_flags if f.cluster_idx == c_num]
        if not matching_flags or matching_flags[0].n_candidates == 0:
            return details

        flag = matching_flags[0]
        # Get matching candidate hits for this cluster
        matching_hits = [
            h for h in essential_hits
            if h.cluster_idx == c_num and h.gene_id in flag.gene_ids
        ]

        if not matching_hits:
            return details

        container = pq("<div class='recruitment-miner-container'>")
        h4 = pq("<h4>")
        h4.text("Target-guided mining: essential-gene paralog(s)")
        container.append(h4)

        summary_p = pq("<p>")
        summary_p.text(f"Prioritization signal: {flag.bonus_signal}")
        container.append(summary_p)

        # Candidate summary lines
        ul = pq("<ul class='recruitment-summary-list'>")
        for h in matching_hits:
            li = pq("<li>")
            gene_safe = html.escape(str(h.gene_id))
            fam_safe = html.escape(str(h.family))
            src_safe = html.escape(str(h.source))
            dup_safe = html.escape(str(h.dup_outside))
            copies_safe = int(h.copies_outside)
            corr_safe = html.escape(str(h.corrupted))
            line_html = f"<b>{gene_safe}</b> — {fam_safe} ({src_safe}), dup outside cluster: {dup_safe} ({copies_safe} copies), corrupted: {corr_safe}"
            li.html(line_html)
            ul.append(li)
        container.append(ul)

        # Collapsible details table
        details_elem = pq("<details class='recruitment-details-block'>")
        summary_elem = pq("<summary>View full alignment details</summary>")
        details_elem.append(summary_elem)

        table = pq("<table class='recruitment-table table table-striped table-bordered' style='margin-top:8px;'>")
        thead = pq("<thead><tr><th>Gene ID</th><th>Target Hit</th><th>Source</th><th>Family</th><th>% Ident</th><th>% Cov</th><th>E-value</th><th>Dup Outside</th><th>Corrupted</th><th>Notes</th></tr></thead>")
        table.append(thead)

        tbody = pq("<tbody>")
        for h in matching_hits:
            tr = pq("<tr>")
            tr.append(pq("<td>").text(str(h.gene_id)))
            tr.append(pq("<td>").text(str(h.agi_hit)))
            tr.append(pq("<td>").text(str(h.source)))
            tr.append(pq("<td>").text(str(h.family)))
            tr.append(pq("<td>").text(f"{h.pident:.1f}%"))
            tr.append(pq("<td>").text(f"{h.qcov:.1f}%"))
            tr.append(pq("<td>").text(f"{h.evalue:.2e}"))
            tr.append(pq("<td>").text(f"{h.dup_outside} ({h.copies_outside} copies)"))
            tr.append(pq("<td>").text(str(h.corrupted)))
            tr.append(pq("<td>").text(str(h.corruption_reason)))
            tbody.append(tr)

        table.append(tbody)
        details_elem.append(table)
        container.append(details_elem)

        if details is None:
            details = pq("<div>")
        details.append(container)
        return details

    except Exception as e:
        logging.warning("recruitment_miner: HTML generation failed: %s", e)
        return details
