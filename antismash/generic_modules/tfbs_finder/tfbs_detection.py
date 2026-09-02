
#
# Copyright (C) 2025 Hannah E. Augustijn 
# Wageningen University & Research & Leiden University
# Department: Department of Bioinformatics & Institute of Biology Leiden
#
# Copyright (C) 2025 Elena Del Pupo
# Wageningen University & Research
# Bioinformatics Group
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
TFBS detection (per-BGC) using PWMs with MOODS.

This implementation:
- iterates over each BGC (cluster) on the record,
- builds ±range promoter windows for CDS that overlap the BGC,
- clips each window to the BGC span,
- merges overlapping windows within that BGC,
- scans each merged interval exactly once,
- aggregates hits across all BGCs (the HTML/output module maps hits to clusters).

Public entry point: run_tfbs_finder(record, pvalue, start_overlap, matrix_path=PWM_PATH)
"""

from __future__ import annotations
from dataclasses import dataclass
from enum import IntEnum, auto
from typing import Any, Dict, List, Optional, Tuple
import bisect
import os
import json
import logging
import multiprocessing
import tempfile
import time

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation

# MOODS is optional (C extension; may be absent if its build failed).
# tfbs_finder degrades gracefully when missing - see run_tfbs_finder.
try:
    from MOODS import tools as moods_tools
    from MOODS import scan, parsers
    _HAS_MOODS = True
except Exception:  # pragma: no cover - depends on env
    moods_tools = None
    scan = None
    parsers = None
    _HAS_MOODS = False

_MOODS_WARNED = False

from antismash import utils

# --------------------------------------------------------------------------
# Constants / cache / multiprocessing global
# --------------------------------------------------------------------------

PWM_PATH = utils.get_full_path(__file__, os.path.join("data", "Athaliana_motifs.manualfromexcluded.json"))
_MATRIX_CACHE: Dict[str, List["Matrix"]] = {}   # cache parsed matrices per file path
_POOL_MATRICES: Optional[List["Matrix"]] = None  # inherited copy-on-write by worker processes after fork


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------

class Confidence(IntEnum):
    WEAK = auto()
    MEDIUM = auto()
    STRONG = auto()

    def __str__(self) -> str:
        return self.name.lower()


@dataclass
class Matrix:
    name: str
    pwm: List[List[float]]           # 4×N, PFM or log-odds
    max_score: float
    min_score: float
    description: str
    species: str
    link: str
    consensus: str
    _threshold: float = -1.0
    is_log_odds: bool = False        # if True, pwm already is log-odds

    @property
    def score_threshold(self) -> float:
        if self._threshold < 0:
            self._threshold = (self.min_score + self.max_score) / 2
        return self._threshold

    def get_score_confidence(self, score: float) -> Confidence:
        if score <= self.min_score:
            return Confidence.WEAK
        if score >= self.score_threshold:
            return Confidence.STRONG
        return Confidence.MEDIUM

    def to_json(self) -> Dict[str, Any]:
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    @staticmethod
    def from_json(name: str, data: Dict[str, Any]) -> "Matrix":
        return Matrix(
            name=name,
            pwm=data["pwm"],
            max_score=data.get("max_score", 0.0),
            min_score=data.get("min_score", 0.0),
            description=data.get("description", ""),
            species=data.get("species", ""),
            link=data.get("link", ""),
            consensus=data.get("consensus", ""),
            is_log_odds=data.get("is_log_odds", False),
        )


@dataclass
class TFBSHit:
    name: str
    start: int                  # absolute genomic coord on record (0-based)
    species: str
    link: str
    description: str
    consensus: str
    confidence: Confidence
    strand: int                 # +1 / -1
    score: float
    max_score: float

    def to_json(self) -> Dict[str, Any]:
        data = dict(vars(self))
        data["confidence"] = str(data["confidence"])
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "TFBSHit":
        d = dict(data)
        d["confidence"] = Confidence[d["confidence"].upper()]
        return TFBSHit(**d)


class TFBSFinderResults:
    schema_version = 1

    def __init__(self, record_id: str, pvalue: float, start_overlap: int,
                 hits_by_record: Dict[str, List[TFBSHit]]) -> None:
        self.record_id = record_id
        self.pvalue = pvalue
        self.start_overlap = start_overlap
        self.hits_by_record = hits_by_record

    def to_json(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "pvalue": self.pvalue,
            "start_overlap": self.start_overlap,
            "hits_by_record": {
                k: [hit.to_json() for hit in v] for k, v in self.hits_by_record.items()
            }
        }

    def get_hits_for_record(self, record_id: str,
                            confidence: Optional[Confidence] = None,
                            allow_better: bool = False) -> List[TFBSHit]:
        hits = self.hits_by_record.get(record_id, [])
        if confidence is None:
            return hits
        if allow_better:
            return [h for h in hits if h.confidence >= confidence]
        return [h for h in hits if h.confidence == confidence]

    def format_html(self) -> str:
        out = [f"<h3>TFBS Finder Results</h3>"]
        if not self.hits_by_record:
            out.append("<p>No transcription factor binding sites detected.</p>")
            return "\n".join(out)
        out += [
            "<table class='table table-sm'>",
            "<thead><tr><th>Motif</th><th>Start</th><th>Strand</th><th>Score</th>"
            "<th>Confidence</th><th>Species</th></tr></thead>",
            "<tbody>",
        ]
        for _, hits in self.hits_by_record.items():
            for h in hits:
                strand = "+" if h.strand == 1 else "−"
                out.append(
                    f"<tr><td>{h.name}</td><td>{h.start}</td><td>{strand}</td>"
                    f"<td>{h.score:.2f}/{h.max_score:.2f}</td>"
                    f"<td>{str(h.confidence).capitalize()}</td><td>{h.species}</td></tr>"
                )
        out += ["</tbody></table>"]
        return "\n".join(out)

    @staticmethod
    def from_json(previous: Dict[str, Any], record: SeqRecord) -> Optional["TFBSFinderResults"]:
        try:
            if previous.get("schema_version") != TFBSFinderResults.schema_version:
                return None
            if previous.get("record_id") != record.id:
                return None
            pvalue = float(previous["pvalue"])
            start_overlap = int(previous["start_overlap"])
            hits_by_record: Dict[str, List[TFBSHit]] = {}
            for k, hits in previous.get("hits_by_record", {}).items():
                hits_by_record[str(k)] = [TFBSHit.from_json(h) for h in hits]
            return TFBSFinderResults(
                record_id=previous["record_id"],
                pvalue=pvalue,
                start_overlap=start_overlap,
                hits_by_record=hits_by_record,
            )
        except Exception:
            return None


# --------------------------------------------------------------------------
# Helpers: windows, matrices, MOODS
# --------------------------------------------------------------------------

def _cds_tss_and_strand(cds) -> Tuple[Optional[int], Optional[int]]:
    """Return TSS (on forward axis) and strand for a CDS (handles split genes)."""
    loc = getattr(cds, "location", None)
    if loc is None:
        return None, None
    strand = int(getattr(loc, "strand", 1) or 1)
    if isinstance(loc, CompoundLocation) and loc.parts:
        first, last = loc.parts[0], loc.parts[-1]
        tss = int(first.start) if strand == 1 else int(last.end) - 1
    else:
        tss = int(loc.start) if strand == 1 else int(loc.end) - 1
    return tss, strand


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping [a,b] inclusive intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [[intervals[0][0], intervals[0][1]]]
    for s, e in intervals[1:]:
        if s > merged[-1][1] + 1:
            merged.append([s, e])
        else:
            if e > merged[-1][1]:
                merged[-1][1] = e
    return [(int(s), int(e)) for s, e in merged]


def build_upstream_caps(all_cds: List) -> Tuple[List[int], List[int]]:
    """
    Precompute sorted CDS body edges for strand-agnostic promoter capping.

    Returns (body_ends, body_starts): sorted lists of all CDS end coordinates
    (0-based exclusive) and start coordinates (0-based inclusive) on the record,
    regardless of strand. A promoter window may not extend into any neighboring
    CDS body: for a TSS at t, a + strand window is bounded by the largest
    body_end <= t, and a - strand window by the smallest body_start > t.
    """
    body_ends = sorted(int(c.location.end) for c in all_cds)
    body_starts = sorted(int(c.location.start) for c in all_cds)
    return body_ends, body_starts


def get_cds_capped_promoter_window(cds,
                                   caps: Tuple[List[int], List[int]],
                                   upstream_bp: int,
                                   seqlen: int) -> Optional[Tuple[int, int]]:
    """
    Strand-aware promoter window for one CDS, capped against neighboring CDS
    bodies irrespective of their strand (an upstream gene on the opposite strand
    bounds the intergenic region just as well):
      + strand: [max(TSS - upstream_bp, largest body end <= TSS), TSS + 50]
      - strand: [TSS - 50, min(TSS + upstream_bp, smallest body start > TSS - 1)]
    Clipped to contig bounds. Returns inclusive 0-based (a, b), or None.

    Note: TSS is a gene-start proxy from CDS annotation (5' UTR introns are
    ignored); the +/-50 flank around the TSS keeps a sliver of the gene's own
    5' end, matching the original asymmetric window design.
    """
    tss, strand = _cds_tss_and_strand(cds)
    if tss is None or strand is None:
        return None
    body_ends, body_starts = caps
    if strand == 1:
        i = bisect.bisect_right(body_ends, tss)
        cap = body_ends[i - 1] if i > 0 else 0
        a = max(tss - upstream_bp, cap)
        b = tss + 50
    else:
        i = bisect.bisect_left(body_starts, tss + 1)
        cap = body_starts[i] - 1 if i < len(body_starts) else seqlen - 1
        a = tss - 50
        b = min(tss + upstream_bp, cap)
    a = max(0, a)
    b = min(seqlen - 1, b)
    if b >= a:
        return int(a), int(b)
    return None


def _collect_windows_for_cluster(record: SeqRecord,
                                 cluster_feature,
                                 upstream_bp: int) -> Tuple[List[Tuple[int, int]], int]:
    """
    Build strand-aware promoter windows centered at the TSS proxy (CDS start fallback):
      + strand: [TSS - upstream_bp, TSS + 50]
      - strand: [TSS - 50, TSS + upstream_bp]
    Each window is capped at neighboring CDS body edges (any strand), i.e. it
    never extends into an upstream neighbor's coding sequence.
    Then clip each window to the cluster span and contig bounds.

    Note: TSS is a gene-start proxy based on CDS annotation (5'-UTR introns are ignored).

    Returns (intervals_inclusive, cds_count_included).
    """
    seqlen = len(record.seq)
    cstart = int(cluster_feature.location.start)
    cend   = int(cluster_feature.location.end) - 1  # inclusive
    cds_count = 0
    raw: List[Tuple[int, int]] = []

    all_cds = list(utils.get_cds_features(record))
    caps = build_upstream_caps(all_cds)

    for cds in all_cds:
        cds_start = int(cds.location.start)
        cds_end   = int(cds.location.end) - 1  # inclusive
        # only consider CDS that overlap this cluster span
        if cds_end < cstart or cds_start > cend:
            continue

        w = get_cds_capped_promoter_window(cds, caps, upstream_bp, seqlen)
        if w is None:
            continue

        # clip to cluster
        a = max(w[0], cstart)
        b = min(w[1], cend)

        if b >= a:
            raw.append((int(a), int(b)))
            cds_count += 1

    return raw, cds_count


def _safe_bg_from_seq(seq: Any) -> List[float]:
    s = str(seq).upper()
    nA = s.count("A")
    nC = s.count("C")
    nG = s.count("G")
    nT = s.count("T")
    total = nA + nC + nG + nT
    if total == 0:
        return [0.25, 0.25, 0.25, 0.25]
    eps = 1e-9
    arr = np.array([nA, nC, nG, nT], dtype=float) + eps
    arr /= arr.sum()
    return arr.tolist()  # A,C,G,T


def _matrix_to_log_odds(matrix: Matrix, background: List[float]) -> List[List[float]]:
    """
    Produce a 4×N log-odds matrix for a motif.
    If matrix.is_log_odds == True, assume matrix.pwm already is log-odds.
    Otherwise, treat matrix.pwm as PFM and convert using MOODS with the *given* background.
    """
    pwm = matrix.pwm if not isinstance(matrix.pwm, np.ndarray) else matrix.pwm.tolist()
    if not pwm or len(pwm) != 4 or any(len(r) != len(pwm[0]) for r in pwm):
        raise ValueError(f"{matrix.name}: PWM must be 4×N")

    # Heuristic: negatives → already log-odds
    if matrix.is_log_odds or any(val < 0 for row in pwm for val in row):
        lod = pwm
    else:
        # Convert PFM -> log-odds via MOODS; needs a temp file
        pfm_str = "\n".join(" ".join(f"{v:.6f}" for v in row) for row in pwm) + "\n"
        tmp = None
        try:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pfm", delete=False)
            tmp.write(pfm_str)
            tmp.flush()
            tmp.close()
            lod = parsers.pfm_to_log_odds(tmp.name, background, 1e-3)
            lod = lod.tolist() if isinstance(lod, np.ndarray) else lod
        finally:
            if tmp is not None:
                try:
                    os.unlink(tmp.name)
                except Exception:
                    pass

    if not lod or len(lod) != 4 or any(len(r) != len(lod[0]) for r in lod):
        raise ValueError(f"{matrix.name}: log-odds must be 4×N")
    return lod


def _load_matrices_cached(json_file: str) -> List[Matrix]:
    mats = _MATRIX_CACHE.get(json_file)
    if mats is not None:
        return mats
    with open(json_file, encoding="utf-8") as fh:
        data = json.load(fh)
    mats = []
    for name, values in data.items():
        try:
            m = Matrix.from_json(name, values)
            if len(m.pwm) != 4 or any(len(row) != len(m.pwm[0]) for row in m.pwm):
                raise ValueError("PWM must be 4×N")
            mats.append(m)
        except Exception as e:
            logging.error("Skipping motif %r due to parse/shape error: %s", name, e)
    _MATRIX_CACHE[json_file] = mats
    logging.debug("Loaded %d matrices from %s", len(mats), json_file)
    return mats


def _absolute_hit(seg_a: int,
                  motif_idx: int,
                  rel_pos: int,
                  strand: int,
                  score: float,
                  fwd_len: int = 0,
                  motif_len: int = 0) -> Tuple[int, int, int, float]:
    """
    Compute absolute coordinate for a local relative hit.
    For strand=+1, absolute pos = seg_a + rel_pos.
    For strand=-1, absolute pos = seg_a + (fwd_len - rel_pos - motif_len).
    """
    if strand == 1:
        abs_pos = seg_a + rel_pos
    else:
        abs_pos = seg_a + (fwd_len - rel_pos - motif_len)
    return (int(motif_idx), int(abs_pos), int(strand), float(score))


def _scan_segment_str_with_pwms(seg_a: int,
                                seg_b: int,
                                seq_str: str,
                                matrices: List[Matrix],
                                pvalue: float) -> List[Tuple[int, int, int, float]]:
    """
    Scan a DNA sequence string representing slice [seg_a:seg_b+1] across PWMs.
    Returns raw hits:
      List[(matrix_idx, absolute_start, strand(+1/-1), score)]
    Uses per-segment background and per-motif MOODS thresholds from p-value.
    """
    fwd_seq = str(seq_str)
    background = _safe_bg_from_seq(fwd_seq)  # A,C,G,T
    hits: List[Tuple[int, int, int, float]] = []
    rc_seq = str(Seq(fwd_seq).reverse_complement())

    # Debug throttle (optional): limit number of motifs via env
    max_pwms = int(os.environ.get("TFBS_MAX_MOTIFS", "0") or "0")
    mats = matrices[:max_pwms] if max_pwms > 0 else matrices

    for idx, m in enumerate(mats):
        try:
            lod = _matrix_to_log_odds(m, background)
            motif_len = len(lod[0])

            thr = moods_tools.threshold_from_p(lod, background, pvalue)
            thresholds = [thr]

            fwd = scan.scan_dna(fwd_seq, [lod], background, thresholds, 7)[0]
            for mh in fwd:
                hits.append(_absolute_hit(seg_a, idx, mh.pos, 1, mh.score))

            rev = scan.scan_dna(rc_seq, [lod], background, thresholds, 7)[0]
            for mh in rev:
                hits.append(_absolute_hit(seg_a, idx, mh.pos, -1, mh.score, len(fwd_seq), motif_len))

        except Exception as e:
            logging.error("MOODS failed for %s: %s", m.name, e)

    return hits


def _scan_segment_with_pwms(record: SeqRecord,
                            seg_a: int,
                            seg_b: int,
                            matrices: List[Matrix],
                            pvalue: float) -> List[Tuple[int, int, int, float]]:
    """
    Scan record.seq[seg_a:seg_b+1] once across all PWMs and return raw hits:
      List[(matrix_idx, absolute_start, strand(+1/-1), score)]
    Uses per-segment background and per-motif MOODS thresholds from p-value.
    """
    seq_str = str(record.seq[seg_a:seg_b+1])
    return _scan_segment_str_with_pwms(seg_a, seg_b, seq_str, matrices, pvalue)


def _scan_interval_worker(args: Tuple[int, int, str, float]) -> List[Tuple[int, int, int, float]]:
    """
    Worker task scanning one interval with MOODS; returns raw hits.
    Reads preloaded matrices from module-level `_POOL_MATRICES` inherited via fork.
    """
    global _POOL_MATRICES
    a, b, seq_str, pvalue = args
    matrices = _POOL_MATRICES if _POOL_MATRICES is not None else []
    return _scan_segment_str_with_pwms(a, b, seq_str, matrices, pvalue)


def _build_worker_tasks(record: SeqRecord,
                        intervals: List[Tuple[int, int]],
                        pvalue: float) -> List[Tuple[int, int, str, float]]:
    """
    Build lightweight per-interval worker args: (a, b, seq_slice_str, pvalue).
    Avoids passing SeqRecord or large parent object graphs to workers.
    """
    return [(int(a), int(b), str(record.seq[a:b+1]), float(pvalue)) for a, b in intervals]


def scan_intervals_parallel(record: SeqRecord,
                            intervals: List[Tuple[int, int]],
                            matrices: List[Matrix],
                            pvalue: float,
                            cpus: int = 1,
                            desc: str = "") -> List[Tuple[int, int, int, float]]:
    """
    Scan a list of (start, end) intervals across PWMs in parallel or sequential mode.

    Returns raw hits in deterministic order:
      List[(matrix_idx, absolute_start, strand(+1/-1), score)]
      sorted by interval start ascending, then motif index / match order.

    Determinism:
      - Sequential path preserves exact interval evaluation order.
      - Parallel path maps intervals with imap preserving input interval order,
        then merges per-interval hit lists.
      - Output intervals are sorted prior to scanning, ensuring identical ordering
        regardless of parallel or sequential execution.

    Fork safety & Memory/IPC efficiency:
      - Slices strings in parent before pool creation; worker args contain only
        (a, b, seq_slice_str, pvalue) (~2 KB per interval, no SeqRecord).
      - `matrices` stored in module-level `_POOL_MATRICES` before Pool creation,
        inherited copy-on-write across workers via multiprocessing fork context.
      - If Pool creation/mapping fails (e.g. non-fork platform or environment restriction),
        gracefully logs a warning and falls back to sequential execution.
    """
    if not intervals or not _HAS_MOODS or not matrices:
        return []

    global _POOL_MATRICES
    _POOL_MATRICES = matrices

    # Sort intervals for deterministic processing
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    n_intervals = len(sorted_intervals)

    # Use sequential loop if cpus <= 1, small number of intervals, or MOODS unavailable
    if cpus <= 1 or n_intervals <= 3:
        all_hits: List[Tuple[int, int, int, float]] = []
        prog_interval = max(1, n_intervals // 10)
        prog_interval = min(prog_interval, 20)
        scanned_bp = 0
        t_seq_start = time.time()
        for j, (a, b) in enumerate(sorted_intervals, 1):
            seg_hits = _scan_segment_with_pwms(record, a, b, matrices, pvalue)
            all_hits.extend(seg_hits)
            scanned_bp += (b - a + 1)
            if j % prog_interval == 0 or j == n_intervals:
                elapsed = time.time() - t_seq_start
                prefix = f"TFBS: {desc} " if desc else "TFBS: "
                logging.warning("%sscanned %d/%d intervals (%.1f%%), bp=%d, hits=%d, elapsed=%.1fs",
                                prefix, j, n_intervals, 100.0 * j / n_intervals,
                                scanned_bp, len(all_hits), elapsed)
        return all_hits

    # Parallel path: slice small string payloads in parent to avoid shipping SeqRecord / large graphs
    t_start = time.time()
    n_workers = min(int(cpus), n_intervals)
    chunksize = max(1, n_intervals // (n_workers * 8))

    tasks = _build_worker_tasks(record, sorted_intervals, pvalue)

    try:
        # Fork-safety: matrices are loaded in the PARENT before Pool creation (inherited copy-on-write)
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(processes=n_workers) as pool:
            results = list(pool.imap(_scan_interval_worker, tasks, chunksize=chunksize))

        all_hits = []
        for seg_hits in results:
            all_hits.extend(seg_hits)

        elapsed = time.time() - t_start
        desc_str = f" ({desc})" if desc else ""
        logging.warning("parallel scan%s: %d intervals on %d workers, %d hits, %.1fs",
                        desc_str, n_intervals, n_workers, len(all_hits), elapsed)
        return all_hits

    except Exception as e:
        logging.warning("TFBS: parallel pool failed%s (%s); falling back to sequential scan",
                        f" for {desc}" if desc else "", e)
        all_hits = []
        for a, b in sorted_intervals:
            seg_hits = _scan_segment_with_pwms(record, a, b, matrices, pvalue)
            all_hits.extend(seg_hits)
        return all_hits


def _filter_hits_to_objects(matrices: List[Matrix],
                            raw_hits: List[Tuple[int, int, int, float]]) -> List[TFBSHit]:
    out: List[TFBSHit] = []
    for mat_idx, start, strand, score in raw_hits:
        m = matrices[mat_idx]
        conf = m.get_score_confidence(score)
        out.append(
            TFBSHit(
                name=m.name,
                start=int(start),
                species=m.species,
                link=m.link,
                description=m.description,
                consensus=m.consensus,
                confidence=conf,
                strand=int(strand),
                score=float(score),
                max_score=float(m.max_score),
            )
        )
    return out


# --------------------------------------------------------------------------
# Public entry point (per-BGC scanning only)
# --------------------------------------------------------------------------

def run_tfbs_finder(record: SeqRecord,
                    pvalue: float,
                    start_overlap: int,
                    matrix_path: str = PWM_PATH,
                    cpus: int = 1) -> TFBSFinderResults:
    """
    Run TFBS scan **per BGC** on this record:
      - for each cluster on the record, build ±start_overlap windows around CDS TSS,
        clipped to the cluster span;
      - merge windows inside that cluster and scan each merged interval once;
      - aggregate/deduplicate hits across clusters.

    Returns TFBSFinderResults with hits_by_record = { record.id: [TFBSHit, ...] }.
    """
    logging.info("TFBS: %s starting (per-BGC mode)", record.id)

    global _MOODS_WARNED
    if not _HAS_MOODS:
        if not _MOODS_WARNED:
            logging.warning("TFBS: MOODS is unavailable (C extension missing); TFBS results will be empty. "
                            "Install MOODS-python via 'pip install MOODS-python' or extra 'pip install -e .[tfbs]'")
            _MOODS_WARNED = True
        return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: []})
    matrices = _load_matrices_cached(matrix_path)
    if not matrices:
        logging.warning("TFBS: no matrices loaded (%s)", matrix_path)
        return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: []})

    clusters = utils.get_sorted_cluster_features(record)
    if not clusters:
        logging.info("TFBS: %s has no clusters; nothing to scan", record.id)
        return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: []})

    all_raw_hits: List[Tuple[int, int, int, float]] = []
    total_bp = 0
    total_int = 0
    t_scan_start = time.time()

    total_raw_hits = 0
    cluster_hit_counts: List[Tuple[int, int]] = []
    for c in clusters:
        cidx = utils.get_cluster_number(c)
        cstart = int(c.location.start)
        cend   = int(c.location.end) - 1  # inclusive

        raw_windows, cds_count = _collect_windows_for_cluster(record, c, start_overlap)
        if not raw_windows:
            logging.warning("TFBS: cluster #%d %d-%d: no CDS windows; skip",
                            cidx, cstart, cend)
            continue

        merged = _merge_intervals(raw_windows)
        bp = sum(b - a + 1 for a, b in merged)
        total_bp += bp
        total_int += len(merged)
        logging.warning("TFBS: cluster #%d %d-%d: CDS windows=%d; merged intervals=%d; bp=%d",
                        cidx, cstart, cend, cds_count, len(merged), bp)

        # Optional debug throttle: limit intervals per cluster
        max_intervals = int(os.environ.get("TFBS_MAX_INTERVALS", "0") or "0")
        intervals = merged[:max_intervals] if max_intervals > 0 else merged

        cluster_hits = scan_intervals_parallel(
            record=record,
            intervals=intervals,
            matrices=matrices,
            pvalue=pvalue,
            cpus=cpus,
            desc=f"cluster #{cidx}",
        )
        total_raw_hits += len(cluster_hits)
        cluster_hit_counts.append((cidx, len(cluster_hits)))
        all_raw_hits.extend(cluster_hits)

    # Diagnostic check for repetitive-promoter hit inflation
    if cluster_hit_counts:
        counts_only = sorted(cnt for _, cnt in cluster_hit_counts)
        n_c = len(counts_only)
        if n_c % 2 == 1:
            median_hits = float(counts_only[n_c // 2])
        else:
            median_hits = (counts_only[n_c // 2 - 1] + counts_only[n_c // 2]) / 2.0
        if median_hits > 50:
            for cidx, c_cnt in cluster_hit_counts:
                if c_cnt > 8.0 * median_hits:
                    logging.warning(
                        "TFBS: cluster #%d hit count %d is %.1fx the record median (%d) — possible repetitive-promoter inflation",
                        cidx,
                        c_cnt,
                        c_cnt / median_hits,
                        int(round(median_hits)),
                    )

    t_scan_elapsed = time.time() - t_scan_start
    logging.warning("⏱ TFBS stage 'cluster_scan' took %.2fs (%d intervals, %d bp, %d raw hits)",
                    t_scan_elapsed, total_int, total_bp, total_raw_hits)

    # Deduplicate across clusters (same motif, start, strand, score)
    if all_raw_hits:
        as_set = {}
        for mi, st, sd, sc in all_raw_hits:
            as_set[(mi, st, sd, round(sc, 4))] = (mi, st, sd, sc)
        all_raw_hits = list(as_set.values())

    logging.warning("🧪 Filtering hits")
    t_filter_start = time.time()
    hits = _filter_hits_to_objects(matrices, all_raw_hits)
    t_filter_elapsed = time.time() - t_filter_start
    logging.warning("⏱ TFBS stage 'hit_filtering' took %.2fs (%d hits)",
                    t_filter_elapsed, len(hits))

    # Attach to record.annotations for downstream HTML (best effort)
    try:
        record.annotations.setdefault("tfbs_finder", {})[record.id] = [h.to_json() for h in hits]
    except Exception as e:
        logging.debug("TFBS: could not attach hits to record annotations: %s", e)

    logging.info("TFBS: %s finished with %d hit record(s)", record.id, len(hits))
    return TFBSFinderResults(record.id, pvalue, start_overlap, {record.id: hits})