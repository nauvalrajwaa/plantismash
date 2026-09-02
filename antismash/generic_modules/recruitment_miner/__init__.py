"""
recruitment_miner — Target-guided mining & primary metabolism paralog detection.

Identifies essential/primary-metabolism gene paralogs inside candidate BGCs
outside the biosynthetic quorum, plus corruption/pseudogenization tags,
as a prioritization signal.
"""

from argparse import Namespace
import logging
import os
import sys
from typing import Any, List, Optional

from Bio.SeqRecord import SeqRecord
from antismash import utils

from .detect import (
    DEFAULT_EVALUE,
    DEFAULT_PIDENT,
    DEFAULT_QCOV,
    EssentialHit,
    ClusterRecruitmentFlag,
    analyze_record_recruitment,
    locate_diamond_binary,
)
from .output import (
    write_recruitment_tsvs,
    generate_details_div,
    _get_output_dir,
)

NAME = "recruitment_miner"
SHORT_DESCRIPTION = "Detects essential/primary-metabolism gene paralogs inside BGCs (target-guided prioritization)"

# Default reference data path
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

_required_binaries = [
    ("diamond", False),
]


def _enabled(options: Any) -> bool:
    return bool(getattr(options, "recruitment_miner", False))


def _resolve_db_dir(options: Any) -> str:
    db_dir = getattr(options, "recruitment_db", None)
    if db_dir and isinstance(db_dir, str) and db_dir.strip():
        return os.path.abspath(db_dir.strip())
    return DATA_DIR


def check_prereqs(options: Any) -> List[str]:
    """Check prerequisites when module is enabled."""
    if not _enabled(options):
        return []

    failures = []
    diamond_bin = locate_diamond_binary()
    if not diamond_bin:
        failures.append("Failed to locate DIAMOND binary required for recruitment_miner")

    db_dir = _resolve_db_dir(options)
    dmnd_file = os.path.join(db_dir, "essential_proteins.dmnd")
    faa_file = os.path.join(db_dir, "essential_proteins.faa")
    if not (os.path.isfile(dmnd_file) or os.path.isfile(faa_file)):
        logging.warning(
            "recruitment_miner: Reference DB not found at %s. Run bash_scripts/build_recruitment_db.py to build.",
            db_dir,
        )

    return failures


def check_options(options: Any) -> List[str]:
    """Validate options."""
    if not _enabled(options):
        return []
    return []


def run_recruitment_miner_for_record(record: SeqRecord, options: Any) -> None:
    """
    Main entry point called per SeqRecord.
    Gated by --recruitment-miner; wraps execution in try/except with logging.
    """
    if not _enabled(options):
        return

    # Ensure extrarecord initialized
    if not hasattr(options, "extrarecord") or options.extrarecord is None:
        options.extrarecord = {}
    options.extrarecord.setdefault(record.id, Namespace())
    ns = options.extrarecord[record.id]
    ns.extradata = getattr(ns, "extradata", {})

    db_dir = _resolve_db_dir(options)
    try:
        utils.log_status(
            "➡️ Running recruitment miner for contig #%d" % getattr(options, "record_idx", 0)
        )
    except Exception:
        pass

    essential_hits: List[EssentialHit] = []
    cluster_flags: List[ClusterRecruitmentFlag] = []

    try:
        essential_hits, cluster_flags = analyze_record_recruitment(
            record=record,
            options=options,
            db_dir=db_dir,
            evalue_cutoff=DEFAULT_EVALUE,
            pident_cutoff=DEFAULT_PIDENT,
            qcov_cutoff=DEFAULT_QCOV,
        )
        logging.info(
            "recruitment_miner: record %s finished with %d essential hit(s), %d flagged cluster(s)",
            record.id,
            len(essential_hits),
            len(cluster_flags),
        )
    except Exception:
        logging.exception("💥 recruitment_miner crashed on %s", record.id)
        essential_hits = []
        cluster_flags = []

    # Stash in options.extrarecord for HTML rendering and run finalization
    ns.extradata["RecruitmentEssentialHits"] = essential_hits
    ns.extradata["RecruitmentClusterFlags"] = cluster_flags
    logging.info("✅ recruitment_miner finished for record: %s", record.id)


def finalize_recruitment_run_outputs(seq_records: List[SeqRecord], options: Any) -> None:
    """
    Run-level finalizer for recruitment_miner.
    Collects EssentialHit and ClusterRecruitmentFlag objects across all records
    and writes outdir/recruitment_miner/{essential_hits,cluster_flags}.tsv once.
    """
    if not _enabled(options):
        return

    try:
        outdir = _get_output_dir(options)
        accum_hits: List[EssentialHit] = []
        accum_flags: List[ClusterRecruitmentFlag] = []
        extrarecord = getattr(options, "extrarecord", {}) or {}

        for rec in seq_records:
            rec_id = getattr(rec, "id", None)
            if not rec_id or rec_id not in extrarecord:
                continue
            ns = extrarecord[rec_id]
            extradata = getattr(ns, "extradata", {}) or {}
            rec_hits = extradata.get("RecruitmentEssentialHits", [])
            if isinstance(rec_hits, list):
                accum_hits.extend(rec_hits)
            rec_flags = extradata.get("RecruitmentClusterFlags", [])
            if isinstance(rec_flags, list):
                accum_flags.extend(rec_flags)

        write_recruitment_tsvs(accum_hits, accum_flags, outdir)
        logging.info(
            "recruitment_miner: finalized TSVs across %d record(s) (%d hits, %d flags)",
            len(seq_records),
            len(accum_hits),
            len(accum_flags),
        )
    except Exception:
        logging.exception("💥 recruitment_miner: failed to finalize run outputs")
