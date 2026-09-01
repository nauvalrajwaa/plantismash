#!/usr/bin/env python3
"""
Unit tests for recruitment_miner (Feature B0 + B1).

Tests pure Python logic without requiring DIAMOND execution:
  (a) Corruption detection (synthetic SeqFeatures: /pseudo, internal stop, truncation vs clean)
  (b) Quorum-exclusion logic (sec_met-labeled vs unlabeled)
  (c) Duplication selection from canned DIAMOND tabular strings (parse function tested directly)
  (d) TSV writing (header-only + populated)
  (e) Candidate selection on a synthetic record with 1 cluster + 2 genes
  (f) HTML generator integration
"""

import os
import shutil
import sys
import tempfile
import unittest
from argparse import Namespace
from unittest.mock import patch

# Ensure repo is in sys.path
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from antismash.generic_modules.recruitment_miner.detect import (
    AlignmentMatch,
    ClusterRecruitmentFlag,
    EssentialHit,
    analyze_record_recruitment,
    check_corruption,
    get_cluster_quorum_gene_ids,
    load_reference_metadata,
    parse_diamond_tabular,
)
from antismash.generic_modules.recruitment_miner.output import (
    CLUSTER_FLAGS_HEADER,
    ESSENTIAL_HITS_HEADER,
    generate_details_div,
    write_recruitment_tsvs,
)
from antismash.generic_modules.recruitment_miner import (
    check_prereqs,
    check_options,
    run_recruitment_miner_for_record,
)


class TestRecruitmentMiner(unittest.TestCase):

    def test_corruption_detection(self):
        """Test corruption and pseudogenization detection rules."""
        # 1. Clean feature
        clean_feat = SeqFeature(
            FeatureLocation(100, 400, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["GENE_CLEAN"], "translation": ["MKLLTVAAGSTV*"]},
        )
        is_corr, reason = check_corruption(clean_feat)
        self.assertFalse(is_corr)
        self.assertEqual(reason, "clean")

        # 2. Pseudo qualifier
        pseudo_feat = SeqFeature(
            FeatureLocation(100, 400, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["GENE_PSEUDO"], "pseudo": [""]},
        )
        is_corr, reason = check_corruption(pseudo_feat)
        self.assertTrue(is_corr)
        self.assertIn("/pseudo qualifier", reason)

        # 3. Internal stop codon
        stop_feat = SeqFeature(
            FeatureLocation(100, 400, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["GENE_STOP"], "translation": ["MKLL*VAAGSTV*"]},
        )
        is_corr, reason = check_corruption(stop_feat)
        self.assertTrue(is_corr)
        self.assertIn("internal stop codon", reason)

        # 4. Truncation relative to genome copy outside
        trunc_feat = SeqFeature(
            FeatureLocation(100, 400, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["GENE_TRUNC"], "translation": ["MKLLTVAAGSTV*"]},
        )
        db_match = AlignmentMatch(
            qseqid="GENE_TRUNC",
            sseqid="AT1G01010",
            pident=85.0,
            length=150,
            evalue=1e-30,
            qlen=300,
            slen=300,
        )  # qcov = 50%
        outside_match = AlignmentMatch(
            qseqid="GENE_OUTSIDE",
            sseqid="GENE_TRUNC",
            pident=90.0,
            length=280,
            evalue=1e-80,
            qlen=300,
            slen=300,
        )  # qcov = 93.3%
        is_corr, reason = check_corruption(
            trunc_feat, db_match=db_match, outside_matches=[outside_match]
        )
        self.assertTrue(is_corr)
        self.assertIn("truncated alignment", reason)

    def test_quorum_exclusion_logic(self):
        """Test distinguishing quorum (biosynthetic core) genes from recruited non-core genes."""
        record = SeqRecord(Seq("A" * 5000), id="chr1")

        # Cluster covering 1000..4000
        cluster_feat = SeqFeature(
            FeatureLocation(1000, 4000, strand=1),
            type="cluster",
            qualifiers={"product": ["terpene"], "note": ["Cluster number: 1"]},
        )

        # Quorum CDS with sec_met Type
        cds_core = SeqFeature(
            FeatureLocation(1100, 1800, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["CORE_GENE_1"],
                "sec_met": ["Type: terpene", "bitscore: 350.0"],
                "translation": ["MKLLTVAAGSTV*"],
            },
        )

        # Non-core CDS in cluster (recruitment candidate)
        cds_recruited = SeqFeature(
            FeatureLocation(2000, 2800, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["RECRUITED_GENE_1"],
                "translation": ["MASRRLLQTV*"],
            },
        )

        record.features = [cluster_feat, cds_core, cds_recruited]

        quorum_ids = get_cluster_quorum_gene_ids(cluster_feat, record)
        self.assertIn("CORE_GENE_1", quorum_ids)
        self.assertNotIn("RECRUITED_GENE_1", quorum_ids)

    def test_canned_diamond_tabular_parsing(self):
        """Test DIAMOND blastp tabular parser with canned outputs."""
        canned_output = (
            "GENE_1\tAT1G78500.1\t65.4\t240\t1.2e-45\t250\t252\n"
            "GENE_2\tAT3G25800.1\t35.0\t100\t1.0e-05\t300\t310\n"
            "# comment line\n"
            "GENE_3\tAT5G42600.1\t82.1\t400\t0.0\t410\t420\n"
        )
        matches = parse_diamond_tabular(canned_output)
        self.assertEqual(len(matches), 3)

        m0 = matches[0]
        self.assertEqual(m0.qseqid, "GENE_1")
        self.assertEqual(m0.sseqid, "AT1G78500.1")
        self.assertAlmostEqual(m0.pident, 65.4)
        self.assertEqual(m0.length, 240)
        self.assertAlmostEqual(m0.evalue, 1.2e-45)
        self.assertEqual(m0.qlen, 250)
        self.assertEqual(m0.slen, 252)
        self.assertAlmostEqual(m0.qcov, (240 / 250.0) * 100.0)

    def test_tsv_writing(self):
        """Test writing TSV outputs (header-only and populated)."""
        temp_dir = tempfile.mkdtemp(prefix="test_rec_tsv_")
        try:
            # 1. Empty results
            hits_tsv, flags_tsv = write_recruitment_tsvs([], [], temp_dir)
            self.assertTrue(os.path.isfile(hits_tsv))
            self.assertTrue(os.path.isfile(flags_tsv))

            with open(hits_tsv, "r", encoding="utf-8") as fh:
                hits_lines = fh.read().strip().split("\n")
                self.assertEqual(hits_lines[0].split("\t"), ESSENTIAL_HITS_HEADER)
                self.assertEqual(len(hits_lines), 1)

            with open(flags_tsv, "r", encoding="utf-8") as fh:
                flags_lines = fh.read().strip().split("\n")
                self.assertEqual(flags_lines[0].split("\t"), CLUSTER_FLAGS_HEADER)
                self.assertEqual(len(flags_lines), 1)

            # 2. Populated results
            h = EssentialHit(
                record_id="chr1",
                cluster_idx=1,
                product="terpene",
                gene_id="GENE_PARALOG",
                agi_hit="AT1G78500",
                pident=72.5,
                qcov=95.0,
                evalue=1e-50,
                source="CURATED_FAMILY",
                family="OSC",
                dup_outside="yes",
                copies_outside=2,
                corrupted="no",
                corruption_reason="clean",
                quorum_excluded="yes",
            )
            f = ClusterRecruitmentFlag(
                record_id="chr1",
                cluster_idx=1,
                product="terpene",
                n_candidates=1,
                n_corrupted=0,
                gene_ids=["GENE_PARALOG"],
                bonus_signal="1 paralogs (0 corrupted)",
            )

            write_recruitment_tsvs([h], [f], temp_dir)

            with open(hits_tsv, "r", encoding="utf-8") as fh:
                hits_lines = fh.read().strip().split("\n")
                self.assertEqual(len(hits_lines), 2)
                row = hits_lines[1].split("\t")
                self.assertEqual(row[0], "chr1")
                self.assertEqual(row[3], "GENE_PARALOG")
                self.assertEqual(row[4], "AT1G78500")
                self.assertEqual(row[8], "CURATED_FAMILY")
                self.assertEqual(row[9], "OSC")

            with open(flags_tsv, "r", encoding="utf-8") as fh:
                flags_lines = fh.read().strip().split("\n")
                self.assertEqual(len(flags_lines), 2)
                row = flags_lines[1].split("\t")
                self.assertEqual(row[3], "1")
                self.assertEqual(row[6], "1 paralogs (0 corrupted)")

        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def test_candidate_selection_end_to_end_mocked(self):
        """Test candidate selection on a synthetic record with 1 cluster + 2 genes using mocked alignments."""
        record = SeqRecord(Seq("ATGC" * 1000), id="contig_synthetic")

        # Cluster covering 100..800
        cluster_feat = SeqFeature(
            FeatureLocation(100, 800, strand=1),
            type="cluster",
            qualifiers={"product": ["terpene"], "note": ["Cluster number: 1"]},
        )
        # Gene 1 (Core terpene synthase inside cluster)
        cds_core = SeqFeature(
            FeatureLocation(120, 400, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["GENE_CORE"],
                "sec_met": ["Type: terpene", "bitscore: 450.0"],
                "translation": ["MSSTTVVLAR*"],
            },
        )
        # Gene 2 (Recruited primary metabolism enzyme inside cluster)
        cds_rec = SeqFeature(
            FeatureLocation(500, 750, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["GENE_RECRUITED"],
                "translation": ["MKLLTVAAGSTV*"],
            },
        )
        # Gene 3 (Primary metabolism enzyme copy OUTSIDE cluster)
        cds_out = SeqFeature(
            FeatureLocation(1500, 1800, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": ["GENE_PRIMARY_OUTSIDE"],
                "translation": ["MKLLTVAAGSTV*"],
            },
        )

        record.features = [cluster_feat, cds_core, cds_rec, cds_out]

        # Mock DB matches: GENE_RECRUITED hits AT1G78500 (OSC family)
        mock_db_matches = [
            AlignmentMatch(
                qseqid="GENE_RECRUITED",
                sseqid="AT1G78500",
                pident=75.0,
                length=200,
                evalue=1e-40,
                qlen=210,
                slen=220,
            )
        ]
        # Mock self matches: GENE_RECRUITED hits GENE_PRIMARY_OUTSIDE
        mock_self_matches = [
            AlignmentMatch(
                qseqid="GENE_RECRUITED",
                sseqid="GENE_PRIMARY_OUTSIDE",
                pident=92.0,
                length=210,
                evalue=1e-80,
                qlen=210,
                slen=210,
            )
        ]

        temp_db_dir = tempfile.mkdtemp(prefix="test_rec_db_")
        try:
            # Create synthetic metadata TSV
            meta_path = os.path.join(temp_db_dir, "reference_metadata.tsv")
            with open(meta_path, "w", encoding="utf-8") as fh:
                fh.write("agi_locus\tsource\tfamily\tevidence_note\n")
                fh.write("AT1G78500\tCURATED_FAMILY\tOSC\tSterol OSC primary enzyme\n")

            # Fake .dmnd
            with open(os.path.join(temp_db_dir, "essential_proteins.dmnd"), "w") as fh:
                fh.write("fake_dmnd")

            def mock_blastp(query_faa, db_path, threads=1, evalue=1e-10):
                if "record_self" in db_path:
                    return mock_self_matches
                return mock_db_matches

            with patch("antismash.generic_modules.recruitment_miner.detect.run_diamond_blastp", side_effect=mock_blastp):
                hits, flags = analyze_record_recruitment(
                    record=record,
                    options=Namespace(cpus=1),
                    db_dir=temp_db_dir,
                )

                self.assertEqual(len(hits), 1)
                self.assertEqual(hits[0].gene_id, "GENE_RECRUITED")
                self.assertEqual(hits[0].agi_hit, "AT1G78500")
                self.assertEqual(hits[0].family, "OSC")
                self.assertEqual(hits[0].dup_outside, "yes")
                self.assertEqual(hits[0].copies_outside, 1)
                self.assertEqual(hits[0].corrupted, "no")
                self.assertEqual(hits[0].quorum_excluded, "yes")

                self.assertEqual(len(flags), 1)
                self.assertEqual(flags[0].n_candidates, 1)
                self.assertEqual(flags[0].n_corrupted, 0)
                self.assertEqual(flags[0].gene_ids, ["GENE_RECRUITED"])
                self.assertEqual(flags[0].bonus_signal, "1 paralogs (0 corrupted)")
        finally:
            shutil.rmtree(temp_db_dir, ignore_errors=True)

    def test_html_details_rendering(self):
        """Test generate_details_div rendering and escaping."""
        seq_record = SeqRecord(Seq("A" * 1000), id="contig_1")
        cluster_feature = SeqFeature(
            FeatureLocation(100, 800),
            type="cluster",
            qualifiers={"product": ["alkaloid"], "note": ["Cluster number: 2"]},
        )

        h = EssentialHit(
            record_id="contig_1",
            cluster_idx=2,
            product="alkaloid",
            gene_id="GENE_PARALOG_1",
            agi_hit="AT1G78500",
            pident=68.2,
            qcov=88.5,
            evalue=2.4e-40,
            source="CURATED_FAMILY",
            family="SPDS/PMT",
            dup_outside="yes",
            copies_outside=3,
            corrupted="no",
            corruption_reason="clean",
            quorum_excluded="yes",
        )
        f = ClusterRecruitmentFlag(
            record_id="contig_1",
            cluster_idx=2,
            product="alkaloid",
            n_candidates=1,
            n_corrupted=0,
            gene_ids=["GENE_PARALOG_1"],
            bonus_signal="1 paralogs (0 corrupted)",
        )

        options = Namespace()
        options.extrarecord = {
            "contig_1": Namespace(
                extradata={
                    "RecruitmentEssentialHits": [h],
                    "RecruitmentClusterFlags": [f],
                }
            )
        }

        html_div = generate_details_div(cluster_feature, seq_record, options)
        self.assertIsNotNone(html_div)
        rendered = html_div.outerHtml()
        self.assertIn("Target-guided mining: essential-gene paralog(s)", rendered)
        self.assertIn("GENE_PARALOG_1", rendered)
        self.assertIn("SPDS/PMT", rendered)
        self.assertIn("3 copies", rendered)

    def test_module_lifecycle_degrades_gracefully(self):
        """Test check_prereqs, check_options, and execution when disabled or missing DB."""
        # 1. Disabled options
        opt_disabled = Namespace(recruitment_miner=False)
        self.assertEqual(check_prereqs(opt_disabled), [])
        self.assertEqual(check_options(opt_disabled), [])

        record = SeqRecord(Seq("ATGC" * 100), id="rec_1")
        # Run on disabled should be no-op without modifying extrarecord
        run_recruitment_miner_for_record(record, opt_disabled)

        # 2. Enabled options with empty/missing DB dir
        temp_empty_dir = tempfile.mkdtemp(prefix="test_rec_empty_")
        try:
            opt_enabled = Namespace(
                recruitment_miner=True,
                recruitment_db=temp_empty_dir,
                full_outputfolder_path=temp_empty_dir,
                record_idx=1,
                cpus=1,
            )
            # check_options passes
            self.assertEqual(check_options(opt_enabled), [])

            # run_recruitment_miner_for_record handles missing DB without crashing
            run_recruitment_miner_for_record(record, opt_enabled)

            # Files still created with headers
            self.assertTrue(
                os.path.isfile(os.path.join(temp_empty_dir, "recruitment_miner", "essential_hits.tsv"))
            )
            self.assertTrue(
                os.path.isfile(os.path.join(temp_empty_dir, "recruitment_miner", "cluster_flags.tsv"))
            )
        finally:
            shutil.rmtree(temp_empty_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
