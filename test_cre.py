#!/usr/bin/env python3
"""
Unit tests for Feature A0/A1 (CRE Motif Enrichment & Scanner Hardening).
Tests pure Python statistics, Benjamini-Hochberg FDR, family aggregation math,
and promoter neighbor capping without requiring the optional MOODS C extension.
"""

import math
import sys
import os
import tempfile
from argparse import Namespace

# Ensure repo is importable
sys.path.insert(0, ".")

from antismash.generic_modules.tfbs_finder.enrichment import (
    log_poisson_pmf,
    poisson_upper_tail,
    benjamini_hochberg_fdr,
    sample_background_windows,
    write_cre_enrichment_tsv,
    CREEnrichmentResult,
    CRE_TSV_HEADER,
)
from antismash.generic_modules.tfbs_finder import finalize_tfbs_run_outputs
from antismash.generic_modules.tfbs_finder.tfbs_detection import (
    _merge_intervals,
    _cds_tss_and_strand,
    _collect_windows_for_cluster,
    _build_worker_tasks,
    _absolute_hit,
    build_upstream_caps,
    get_cds_capped_promoter_window,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


def test_poisson_tail():
    # 1. Compare against known analytical values
    # For lambda = 2.0, P(X >= 5) = 1 - sum_{k=0..4} exp(-2)*2^k / k!
    # P(0) = e^-2, P(1)=2e^-2, P(2)=2e^-2, P(3)=4/3 e^-2, P(4)=2/3 e^-2
    # sum(0..4) = (1 + 2 + 2 + 4/3 + 2/3) * e^-2 = 7 * e^-2 = 0.9473469826562897
    # P(X >= 5) = 1 - 7 * e^-2 = 0.0526530173437103
    p_known_1 = 1.0 - 7.0 * math.exp(-2.0)
    p_calc_1 = poisson_upper_tail(5, 2.0)
    assert abs(p_calc_1 - p_known_1) < 1e-9, f"Poisson tail mismatch: {p_calc_1} vs {p_known_1}"

    # For lambda = 0.5, P(X >= 2) = 1 - (P(0) + P(1)) = 1 - (e^-0.5 + 0.5*e^-0.5) = 1 - 1.5 * e^-0.5
    # 1 - 1.5 * exp(-0.5) = 0.09020401043104986
    p_known_2 = 1.0 - 1.5 * math.exp(-0.5)
    p_calc_2 = poisson_upper_tail(2, 0.5)
    assert abs(p_calc_2 - p_known_2) < 1e-9, f"Poisson tail mismatch: {p_calc_2} vs {p_known_2}"

    # Edge cases
    assert poisson_upper_tail(0, 5.0) == 1.0
    assert poisson_upper_tail(10, 0.0) == 0.0

    # k > lambda: direct upper tail must not collapse to exactly 0.0
    p_ext = poisson_upper_tail(30, 5.0)
    assert 0.0 < p_ext < 1e-13, f"extreme tail collapsed: {p_ext}"
    assert poisson_upper_tail(25, 5.0) > p_ext, "tail must be monotone in k"
    # k == lambda boundary uses the subtraction branch, still in (0, 1)
    p_mode = poisson_upper_tail(5, 5.0)
    assert 0.0 < p_mode < 0.7, f"k == lambda misbehaved: {p_mode}"

    print("test_poisson_tail ok")


def test_bh_fdr():
    # Toy p-value vector
    p_vals = [0.01, 0.04, 0.03, 0.20]
    # Sorted: p[0]=0.01 (rank 1), p[2]=0.03 (rank 2), p[1]=0.04 (rank 3), p[3]=0.20 (rank 4)
    # raw q:
    # rank 4 (0.20): 0.20 * 4 / 4 = 0.20
    # rank 3 (0.04): 0.04 * 4 / 3 = 0.05333333333333334
    # rank 2 (0.03): 0.03 * 4 / 2 = 0.06 -> monotonic with rank 3 min(0.06, 0.05333) = 0.05333333333333334
    # rank 1 (0.01): 0.01 * 4 / 1 = 0.04 -> monotonic min(0.04, 0.05333) = 0.04
    q_vals = benjamini_hochberg_fdr(p_vals)

    assert abs(q_vals[0] - 0.04) < 1e-9
    assert abs(q_vals[1] - 0.05333333333333334) < 1e-9
    assert abs(q_vals[2] - 0.05333333333333334) < 1e-9
    assert abs(q_vals[3] - 0.20) < 1e-9

    assert benjamini_hochberg_fdr([]) == []
    print("test_bh_fdr ok")


def test_family_aggregation_and_math():
    # Synthetic family with 3 motifs
    # Background rates (hits/kb): m1=0.5, m2=1.0, m3=0.5 -> family rate = 2.0 hits/kb
    bg_rates = {"m1": 0.5, "m2": 1.0, "m3": 0.5}
    family_motifs = ["m1", "m2", "m3"]

    # Cluster scanned promoter kb = 1.5 kb
    promoter_kb = 1.5
    exp_hits = sum(bg_rates[m] for m in family_motifs) * promoter_kb
    assert abs(exp_hits - 3.0) < 1e-9

    # Observed hits = 7
    obs_hits = 7
    fold_enrichment = obs_hits / exp_hits
    assert abs(fold_enrichment - (7.0 / 3.0)) < 1e-9

    p_val = poisson_upper_tail(obs_hits, exp_hits)
    assert p_val < 0.05

    # Coherence: 3 out of 4 quorum genes have >= 1 hit
    n_quorum = 4
    n_with_hit = 3
    coherence = n_with_hit / float(n_quorum)
    assert abs(coherence - 0.75) < 1e-9

    res = CREEnrichmentResult(
        record="chr1",
        cluster="1",
        product="terpene",
        n_quorum_genes=n_quorum,
        promoter_kb_scanned=promoter_kb,
        tf_family="MYB",
        n_motifs=len(family_motifs),
        hits_obs=obs_hits,
        hits_exp=exp_hits,
        fold_enrichment=fold_enrichment,
        p_value=p_val,
        q_value=p_val,
        coherence_fraction=coherence,
        top_motifs="m1;m2",
        max_confidence="Strong",
    )
    tsv_row = res.to_tsv_row()
    assert len(tsv_row) == len(CRE_TSV_HEADER)
    assert tsv_row[0] == "chr1"
    assert tsv_row[5] == "MYB"
    assert tsv_row[7] == "7"

    print("test_family_aggregation_and_math ok")


def test_neighbor_capping():
    # Construct synthetic SeqRecord (length 10,000 bp)
    seq = Seq("N" * 10000)
    record = SeqRecord(seq, id="test_rec")

    # Correct semantics: a promoter window may NOT extend into ANY neighboring
    # CDS body (strand-agnostic). loc.end is 0-based EXCLUSIVE, so a + strand
    # window is capped at the upstream neighbor's body END (= first intergenic
    # base), and a - strand window at the downstream neighbor's body START - 1.

    # 1. Plus strand CDS 1: body [1000, 2000) strand +1 (TSS = 1000)
    # 2. Plus strand CDS 2: body [2500, 3500) strand +1 (TSS = 2500)
    # CDS 1: no upstream body -> (0, 1050)
    # CDS 2: capped at CDS 1 body end 2000 -> (2000, 2550)
    #   (uncapped would be (1500, 2550); the old buggy cap would be (1000, 2550))
    f_plus1 = SeqFeature(FeatureLocation(1000, 2000, strand=1), type="CDS")
    f_plus2 = SeqFeature(FeatureLocation(2500, 3500, strand=1), type="CDS")

    # Minus strand CDS 1: body [5000, 6000) strand -1 (TSS = 5999)
    # Minus strand CDS 2: body [7000, 8000) strand -1 (TSS = 7999)
    # CDS 1 (upstream_bp=1500): capped before CDS 2 body start 7000 -> (5949, 6999)
    #   (uncapped would be (5949, 7499) — reaching into CDS 2's body)
    # CDS 2 (upstream_bp=1000): no downstream body -> (7949, 8999)
    f_minus1 = SeqFeature(FeatureLocation(5000, 6000, strand=-1), type="CDS")
    f_minus2 = SeqFeature(FeatureLocation(7000, 8000, strand=-1), type="CDS")

    all_cds = [f_plus1, f_plus2, f_minus1, f_minus2]
    caps = build_upstream_caps(all_cds)

    w_p1 = get_cds_capped_promoter_window(f_plus1, caps, 1000, len(seq))
    assert w_p1 == (0, 1050), f"Got {w_p1}"

    w_p2 = get_cds_capped_promoter_window(f_plus2, caps, 1000, len(seq))
    assert w_p2 == (2000, 2550), f"Got {w_p2}"

    w_m1 = get_cds_capped_promoter_window(f_minus1, caps, 1500, len(seq))
    assert w_m1 == (5949, 6999), f"Got {w_m1}"

    w_m2 = get_cds_capped_promoter_window(f_minus2, caps, 1000, len(seq))
    assert w_m2 == (7949, 8999), f"Got {w_m2}"

    # Opposite-strand neighbor also bounds the window:
    # target + CDS TSS=4000, upstream 1000; a - strand body [3200, 3600)
    # -> capped at its body end 3600 -> (3600, 4050) (uncapped would be (3000, 4050))
    f_tgt = SeqFeature(FeatureLocation(4000, 4600, strand=1), type="CDS")
    f_opp = SeqFeature(FeatureLocation(3200, 3600, strand=-1), type="CDS")
    caps2 = build_upstream_caps([f_tgt, f_opp])
    w_opp = get_cds_capped_promoter_window(f_tgt, caps2, 1000, len(seq))
    assert w_opp == (3600, 4050), f"Got {w_opp}"

    # Adjacent genes (neighbor end == TSS): degenerate window = 51 bp at own 5' end
    f_a = SeqFeature(FeatureLocation(100, 200, strand=1), type="CDS")
    f_b = SeqFeature(FeatureLocation(200, 400, strand=1), type="CDS")
    caps3 = build_upstream_caps([f_a, f_b])
    w_adj = get_cds_capped_promoter_window(f_b, caps3, 1000, len(seq))
    assert w_adj == (200, 250), f"Got {w_adj}"

    # CompoundLocation test
    loc_part1 = FeatureLocation(100, 300, strand=1)
    loc_part2 = FeatureLocation(400, 600, strand=1)
    comp_loc = CompoundLocation([loc_part1, loc_part2])
    f_comp = SeqFeature(comp_loc, type="CDS")
    tss, s = _cds_tss_and_strand(f_comp)
    assert tss == 100 and s == 1

    # Cluster feature test via _collect_windows_for_cluster
    cluster_feat = SeqFeature(FeatureLocation(500, 4000), type="cluster")
    record.features = [f_plus1, f_plus2, cluster_feat]
    raw_windows, count = _collect_windows_for_cluster(record, cluster_feat, 1000)
    assert count == 2
    # Windows clipped to cluster span [500, 3999]
    # f_plus1: (500, 1050); f_plus2: (2000, 2550)
    assert (500, 1050) in raw_windows
    assert (2000, 2550) in raw_windows

    print("test_neighbor_capping ok")


def test_background_window_sampling():
    # 1. Windows under budget -> return all, no cap
    windows = [(100, 500), (1000, 1500), (2000, 2200)]  # lengths: 401, 501, 201 -> total 1103 bp
    sampled, total_bp, sampled_bp = sample_background_windows(windows, 2000)
    assert total_bp == 1103
    assert sampled_bp == 1103
    assert sampled == sorted(windows)

    # 2. Total exactly equals budget -> return all
    sampled, total_bp, sampled_bp = sample_background_windows(windows, 1103)
    assert total_bp == 1103
    assert sampled_bp == 1103
    assert sampled == sorted(windows)

    # 3. Windows over budget -> cap honored, sampled_bp >= budget (first window reaching it)
    # Generate 100 intervals of size 1000 bp (total 100,000 bp)
    large_windows = [(i * 2000, i * 2000 + 999) for i in range(100)]
    budget = 30000  # sample ~30 kb
    sampled, total_bp, sampled_bp = sample_background_windows(large_windows, budget)
    assert total_bp == 100000
    assert sampled_bp >= budget
    assert sampled_bp < 40000
    assert len(sampled) < len(large_windows)
    # Window order must be sorted by coordinates
    assert sampled == sorted(sampled)

    # 4. Deterministic across multiple calls and shuffled input order
    import random
    shuffled_windows = list(large_windows)
    random.seed(42)
    random.shuffle(shuffled_windows)
    sampled2, total_bp2, sampled_bp2 = sample_background_windows(shuffled_windows, budget)
    assert sampled == sampled2
    assert total_bp == total_bp2
    assert sampled_bp == sampled_bp2

    # 5. Empty windows or zero budget
    assert sample_background_windows([], 1000) == ([], 0, 0)
    assert sample_background_windows(windows, 0) == (sorted(windows), 1103, 1103)

    print("test_background_window_sampling ok")


def test_finalize_tfbs_run_outputs_multi_record():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create synthetic SeqRecords
        rec1 = SeqRecord(Seq("N" * 1000), id="rec1")
        rec2 = SeqRecord(Seq("N" * 1000), id="rec2")
        seq_records = [rec1, rec2]

        # Fabricate CREEnrichmentResults for rec1 and rec2
        r1_c1 = CREEnrichmentResult(
            record="rec1",
            cluster="1",
            product="plants/plant",
            n_quorum_genes=5,
            promoter_kb_scanned=5.0,
            tf_family="WRKY",
            n_motifs=10,
            hits_obs=20,
            hits_exp=10.0,
            fold_enrichment=2.0,
            p_value=0.001,
            q_value=0.01,
            coherence_fraction=0.8,
            top_motifs="AT1G01;AT1G02",
            max_confidence="Strong",
        )
        r1_c2 = CREEnrichmentResult(
            record="rec1",
            cluster="2",
            product="plants/alkaloid",
            n_quorum_genes=4,
            promoter_kb_scanned=4.0,
            tf_family="MYB",
            n_motifs=8,
            hits_obs=15,
            hits_exp=8.0,
            fold_enrichment=1.875,
            p_value=0.005,
            q_value=0.02,
            coherence_fraction=0.75,
            top_motifs="AT2G01",
            max_confidence="Medium",
        )
        r2_c1 = CREEnrichmentResult(
            record="rec2",
            cluster="1",
            product="plants/saccharide",
            n_quorum_genes=6,
            promoter_kb_scanned=6.0,
            tf_family="bHLH",
            n_motifs=12,
            hits_obs=0,
            hits_exp=5.0,
            fold_enrichment=0.0,
            p_value=1.0,
            q_value=1.0,
            coherence_fraction=0.0,
            top_motifs="-",
            max_confidence="-",
        )

        options = Namespace(
            tfbs_detection=True,
            outputfoldername=tmpdir,
            full_outputfolder_path=tmpdir,
            extrarecord={
                "rec1": Namespace(extradata={"CREEnrichmentResults": [r1_c1, r1_c2]}),
                "rec2": Namespace(extradata={"CREEnrichmentResults": [r2_c1]}),
            },
        )

        finalize_tfbs_run_outputs(seq_records, options)

        tsv_path = os.path.join(tmpdir, "tfbs", "cre_enrichment.tsv")
        assert os.path.isfile(tsv_path), f"File not found: {tsv_path}"

        with open(tsv_path, "r", encoding="utf-8") as fh:
            lines = [line.strip().split("\t") for line in fh if line.strip()]

        # Header + 3 data rows from both records
        assert len(lines) == 4, f"Expected 4 lines (header + 3 rows), got {len(lines)}"
        assert lines[0] == CRE_TSV_HEADER

        # Sort order check: (record, cluster, q_value)
        # Row 1: rec1, cluster 1
        assert lines[1][0] == "rec1" and lines[1][1] == "1" and lines[1][5] == "WRKY"
        # Row 2: rec1, cluster 2
        assert lines[2][0] == "rec1" and lines[2][1] == "2" and lines[2][5] == "MYB"
        # Row 3: rec2, cluster 1
        assert lines[3][0] == "rec2" and lines[3][1] == "1" and lines[3][5] == "bHLH"
        # Check max_confidence for hits_obs=0 is '-'
        assert lines[3][14] == "-"

    print("test_finalize_tfbs_run_outputs_multi_record ok")


def test_parallel_worker_args_and_absolute_coords():
    # 1. Test _build_worker_tasks slicing on a synthetic record
    seq_str = "ACGTACGTNNNNACGT"
    record = SeqRecord(Seq(seq_str), id="test_rec")
    intervals = [(2, 6), (8, 12)]
    pvalue = 0.001

    tasks = _build_worker_tasks(record, intervals, pvalue)
    assert len(tasks) == 2

    # Verify task 0: (a, b, slice_str, pvalue)
    # interval 2..6 inclusive -> length 5 (indices 2,3,4,5,6)
    expected_slice_0 = str(record.seq[2:6+1])
    assert tasks[0] == (2, 6, expected_slice_0, 0.001)
    assert tasks[0][2] == "GTACG"

    # Verify task 1: (8, 12) inclusive -> length 5 (indices 8,9,10,11,12)
    expected_slice_1 = str(record.seq[8:12+1])
    assert tasks[1] == (8, 12, expected_slice_1, 0.001)
    assert tasks[1][2] == "NNNNA"

    # 2. Assert all worker arg elements are primitive types (no SeqRecord or large graphs)
    for task in tasks:
        assert isinstance(task, tuple)
        assert len(task) == 4
        a, b, s, p = task
        assert type(a) is int
        assert type(b) is int
        assert type(s) is str
        assert type(p) is float
        # Ensure no SeqRecord is referenced in task
        assert not hasattr(task, "seq")
        assert not hasattr(s, "seq")

    # 3. Test _absolute_hit coordinate math for +1 and -1 strands
    seg_a = 1000
    motif_idx = 42
    rel_pos = 15
    score = 8.75
    fwd_len = 200
    motif_len = 10

    # Strand +1: abs_pos = seg_a + rel_pos = 1000 + 15 = 1015
    hit_plus = _absolute_hit(seg_a, motif_idx, rel_pos, 1, score, fwd_len, motif_len)
    assert hit_plus == (42, 1015, 1, 8.75)
    assert type(hit_plus[0]) is int
    assert type(hit_plus[1]) is int
    assert type(hit_plus[2]) is int
    assert type(hit_plus[3]) is float

    # Strand -1: abs_pos = seg_a + (fwd_len - rel_pos - motif_len)
    # = 1000 + (200 - 15 - 10) = 1000 + 175 = 1175
    hit_minus = _absolute_hit(seg_a, motif_idx, rel_pos, -1, score, fwd_len, motif_len)
    assert hit_minus == (42, 1175, -1, 8.75)

    print("test_parallel_worker_args_and_absolute_coords ok")


if __name__ == "__main__":
    test_poisson_tail()
    test_bh_fdr()
    test_family_aggregation_and_math()
    test_neighbor_capping()
    test_background_window_sampling()
    test_finalize_tfbs_run_outputs_multi_record()
    test_parallel_worker_args_and_absolute_coords()
    print("all CRE unit tests passed (7/7)")
