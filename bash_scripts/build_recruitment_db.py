#!/usr/bin/env python3
"""
Reference DB builder for recruitment_miner (B0).

Assembles reference dataset of essential and recruited plant gene families:
  1. EMB/SeedGenes embryo-defective essential loci
  2. OGEE Arabidopsis essential gene annotations
  3. Curated recruited / target-pathway enzyme families (OSC, CPS/KSL, HMGR, SPDS/PMT, SAMS, ACCase-CT, tubulins)

Resolves locus IDs against TAIR10 pep FASTA, deduplicates by locus,
writes essential_proteins.faa and reference_metadata.tsv, and builds essential_proteins.dmnd via diamond makedb.
"""

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from typing import Dict, List, Optional, Set, Tuple

# Default URLs for optional remote fetching
URL_SEEDGENES = os.environ.get(
    "PLANTISMASH_RECRUITMENT_SEEDGENES_URL",
    "http://www.seedgenes.org/Download/SeedGenes_Ath_Mutants.txt",
)
URL_OGEE = os.environ.get(
    "PLANTISMASH_RECRUITMENT_OGEE_URL",
    "https://zenodo.org/records/16926250/files/ogee_athaliana_essential.tsv?download=1",
)
URL_TAIR10_PEP = os.environ.get(
    "PLANTISMASH_RECRUITMENT_TAIR10_URL",
    "https://zenodo.org/records/16926250/files/TAIR10_pep_20101214.fasta?download=1",
)

METADATA_HEADER = ["agi_locus", "source", "family", "evidence_note"]


def parse_tair_pep_fasta(fasta_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse TAIR10 protein FASTA.
    Returns mapping: locus (uppercase, e.g. AT1G01010) -> (protein_id, sequence).
    If multiple splice variants exist, picks .1 variant or longest sequence.
    """
    seqs: Dict[str, List[Tuple[str, str]]] = {}
    current_id = None
    current_seq = []

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    prot_id = current_id.split()[0]
                    locus = prot_id.split(".")[0].upper()
                    seq = "".join(current_seq)
                    seqs.setdefault(locus, []).append((prot_id, seq))
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            prot_id = current_id.split()[0]
            locus = prot_id.split(".")[0].upper()
            seq = "".join(current_seq)
            seqs.setdefault(locus, []).append((prot_id, seq))

    best_seqs: Dict[str, Tuple[str, str]] = {}
    for locus, variants in seqs.items():
        # Prefer .1 splice variant
        dot_one = [v for v in variants if v[0].endswith(".1")]
        if dot_one:
            best_seqs[locus] = dot_one[0]
        else:
            # Pick longest
            best_seqs[locus] = max(variants, key=lambda v: len(v[1]))

    return best_seqs


def parse_emb_file(emb_path: str) -> List[Tuple[str, str, str, str]]:
    """Parse EMB locus list (one AGI locus per line, or TSV)."""
    records = []
    with open(emb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            locus = parts[0].strip().split()[0].upper()
            note = parts[1].strip() if len(parts) > 1 else "SeedGenes/EMB embryo-defective locus"
            if locus.startswith("AT") and len(locus) >= 9:
                records.append((locus, "EMB", "EMB_essential", note))
    return records


def parse_ogee_file(ogee_path: str) -> List[Tuple[str, str, str, str]]:
    """Parse OGEE export TSV. Auto-detects locus and essentiality columns."""
    records = []
    with open(ogee_path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return records

        locus_idx = -1
        ess_idx = -1
        for i, col in enumerate(header):
            col_l = col.lower().strip()
            if col_l in ("locus", "gene", "gene_id", "locus_id", "agi_locus", "symbol"):
                locus_idx = i
            elif col_l in ("essentiality", "essential", "status", "phenotype"):
                ess_idx = i

        if locus_idx == -1:
            locus_idx = 0

        for row in reader:
            if not row or len(row) <= locus_idx:
                continue
            locus = row[locus_idx].strip().upper()
            if "." in locus:
                locus = locus.split(".")[0]
            if not (locus.startswith("AT") and len(locus) >= 9):
                continue

            if ess_idx != -1 and len(row) > ess_idx:
                status = row[ess_idx].strip().lower()
                if "essential" not in status and "ess" not in status:
                    continue

            records.append((locus, "OGEE", "OGEE_essential", "OGEE Arabidopsis essential gene"))
    return records


def parse_families_file(fam_path: str) -> List[Tuple[str, str, str, str]]:
    """
    Parse curated families TSV:
    family<TAB>agi_id[;agi_id...]<TAB>[evidence_note]
    or agi_id<TAB>family
    """
    records = []
    with open(fam_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            col1 = parts[0].strip()
            col2 = parts[1].strip()
            note = parts[2].strip() if len(parts) > 2 else "Curated recruited/target family"

            # Check if col1 is locus or family
            if col1.upper().startswith("AT") and len(col1) >= 9:
                locus = col1.upper()
                family = col2
                records.append((locus, "CURATED_FAMILY", family, note))
            else:
                family = col1
                loci = [loc.strip().upper() for loc in col2.replace(",", ";").split(";") if loc.strip()]
                for locus in loci:
                    if "." in locus:
                        locus = locus.split(".")[0]
                    records.append((locus, "CURATED_FAMILY", family, note))
    return records


def try_fetch_url(url: str, timeout: int = 10) -> Optional[bytes]:
    """Attempt HTTP GET with short timeout. Returns bytes or None."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "plantismash-recruitment-builder/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return response.read()
    except Exception as e:
        logging.warning("Fetch failed for %s: %s", url, e)
        return None


def run_diamond_makedb(fasta_path: str, dmnd_path: str) -> bool:
    """Run `diamond makedb --in fasta_path -d dmnd_path`."""
    diamond_bin = shutil.which("diamond")
    if not diamond_bin:
        # Check standard conda / local paths
        candidates = [
            os.path.expanduser("~/miniforge3/envs/plantismash/bin/diamond"),
            os.path.expanduser("~/miniconda3/envs/plantismash/bin/diamond"),
        ]
        for c in candidates:
            if os.path.isfile(c) and os.access(c, os.X_OK):
                diamond_bin = c
                break

    if not diamond_bin:
        logging.error("DIAMOND executable not found on PATH or standard conda envs.")
        return False

    cmd = [diamond_bin, "makedb", "--in", fasta_path, "-d", dmnd_path]
    logging.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        logging.error("diamond makedb failed (code %d):\n%s", res.returncode, res.stderr)
        return False
    logging.info("DIAMOND database successfully created at %s", dmnd_path)
    return True


def build_recruitment_db(
    tair_pep: str,
    out_dir: str,
    emb_file: Optional[str] = None,
    ogee_file: Optional[str] = None,
    families_file: Optional[str] = None,
    fetch_missing: bool = False,
) -> bool:
    """Main build workflow."""
    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse TAIR10 pep
    if not os.path.isfile(tair_pep):
        if fetch_missing:
            logging.info("TAIR10 pep not found locally at %s, attempting fetch from %s", tair_pep, URL_TAIR10_PEP)
            data = try_fetch_url(URL_TAIR10_PEP)
            if data:
                with open(tair_pep, "wb") as fh:
                    fh.write(data)
            else:
                logging.warning("Could not download TAIR10 pep. Exiting.")
                return False
        else:
            logging.error("TAIR10 pep FASTA file not found: %s", tair_pep)
            return False

    logging.info("Parsing TAIR10 peptide sequences from %s...", tair_pep)
    tair_seqs = parse_tair_pep_fasta(tair_pep)
    logging.info("Parsed %d unique AGI loci from TAIR10.", len(tair_seqs))

    # 2. Collect candidate entries
    entries: List[Tuple[str, str, str, str]] = []

    if emb_file and os.path.isfile(emb_file):
        logging.info("Parsing EMB loci from %s...", emb_file)
        entries.extend(parse_emb_file(emb_file))
    elif fetch_missing:
        logging.info("Attempting fetch of EMB loci from %s...", URL_SEEDGENES)
        emb_data = try_fetch_url(URL_SEEDGENES)
        if emb_data:
            with tempfile.NamedTemporaryFile(delete=False, mode="wb") as tf:
                tf.write(emb_data)
                tf_path = tf.name
            entries.extend(parse_emb_file(tf_path))
            os.remove(tf_path)

    if ogee_file and os.path.isfile(ogee_file):
        logging.info("Parsing OGEE entries from %s...", ogee_file)
        entries.extend(parse_ogee_file(ogee_file))
    elif fetch_missing:
        logging.info("Attempting fetch of OGEE entries from %s...", URL_OGEE)
        ogee_data = try_fetch_url(URL_OGEE)
        if ogee_data:
            with tempfile.NamedTemporaryFile(delete=False, mode="wb") as tf:
                tf.write(ogee_data)
                tf_path = tf.name
            entries.extend(parse_ogee_file(tf_path))
            os.remove(tf_path)

    if families_file and os.path.isfile(families_file):
        logging.info("Parsing curated families from %s...", families_file)
        entries.extend(parse_families_file(families_file))

    if not entries:
        logging.warning("No reference entries collected! Check input files.")
        return False

    # 3. Deduplicate by locus and resolve to FASTA sequences
    # Priority: CURATED_FAMILY > EMB > OGEE if multiple sources specify the locus
    source_priority = {"CURATED_FAMILY": 3, "EMB": 2, "OGEE": 1}
    deduped: Dict[str, Tuple[str, str, str]] = {}  # locus -> (source, family, note)

    for locus, source, family, note in entries:
        if locus not in deduped:
            deduped[locus] = (source, family, note)
        else:
            curr_src, curr_fam, curr_note = deduped[locus]
            if source_priority.get(source, 0) > source_priority.get(curr_src, 0):
                deduped[locus] = (source, family, note)
            elif source == curr_src:
                # Merge notes/families if helpful
                if family and family != curr_fam and curr_fam == "EMB_essential":
                    deduped[locus] = (source, family, f"{curr_note}; {note}")

    logging.info("Total unique reference loci: %d", len(deduped))

    # 4. Write essential_proteins.faa and reference_metadata.tsv
    out_faa = os.path.join(out_dir, "essential_proteins.faa")
    out_tsv = os.path.join(out_dir, "reference_metadata.tsv")
    out_dmnd = os.path.join(out_dir, "essential_proteins.dmnd")

    matched_count = 0
    with open(out_faa, "w", encoding="utf-8") as faa_fh, open(out_tsv, "w", encoding="utf-8", newline="") as tsv_fh:
        tsv_writer = csv.writer(tsv_fh, delimiter="\t")
        tsv_writer.writerow(METADATA_HEADER)

        for locus in sorted(deduped.keys()):
            source, family, note = deduped[locus]
            tsv_writer.writerow([locus, source, family, note])

            if locus in tair_seqs:
                prot_id, seq = tair_seqs[locus]
                faa_fh.write(f">{locus} {prot_id} source={source} family={family}\n{seq}\n")
                matched_count += 1
            else:
                logging.debug("Locus %s not found in TAIR10 pep", locus)

    logging.info("Wrote metadata for %d loci to %s", len(deduped), out_tsv)
    logging.info("Wrote %d matching sequences to %s", matched_count, out_faa)

    if matched_count == 0:
        logging.error("No sequences matched between loci and TAIR10 pep.")
        return False

    # 5. Build DIAMOND DB
    ok = run_diamond_makedb(out_faa, out_dmnd)
    return ok


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    default_out = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "antismash", "generic_modules", "recruitment_miner", "data")
    )

    parser = argparse.ArgumentParser(
        description="Build recruitment_miner reference database from EMB/SeedGenes, OGEE, and curated families."
    )
    parser.add_argument("--tair-pep", required=True, help="TAIR10 pep FASTA file path")
    parser.add_argument("--emb-loci", dest="emb_file", default=None, help="EMB / SeedGenes locus file")
    parser.add_argument("--ogee", dest="ogee_file", default=None, help="OGEE essential gene TSV")
    parser.add_argument("--families", dest="families_file", default=None, help="Curated recruited families TSV")
    parser.add_argument("--out-dir", default=default_out, help="Output directory for essential DB")
    parser.add_argument("--fetch", action="store_true", default=False, help="Attempt to fetch missing files from public URLs")

    args = parser.parse_args()

    success = build_recruitment_db(
        tair_pep=args.tair_pep,
        out_dir=args.out_dir,
        emb_file=args.emb_file,
        ogee_file=args.ogee_file,
        families_file=args.families_file,
        fetch_missing=args.fetch,
    )

    if not success:
        logging.warning("Database build did not complete fully. Check inputs or fetch settings.")
        sys.exit(0)
    else:
        logging.info("Recruitment DB build completed successfully.")


if __name__ == "__main__":
    main()
