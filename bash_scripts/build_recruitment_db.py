#!/usr/bin/env python3
"""
Automated Reference DB builder for recruitment_miner (Feature B0/B1).

Assembles reference dataset of essential and recruited plant gene families:
  1. EMB/SeedGenes embryo-defective essential loci
  2. OGEE / Phenotypic Arabidopsis essential gene annotations
  3. Curated recruited / target-pathway enzyme families:
     - OSC (oxidosqualene cyclases / triterpene synthases)
     - CPS / KSL (copalyl diphosphate & kaurene synthases)
     - HMGR (HMG-CoA reductase, mevalonate pathway)
     - SPDS / PMT (spermidine synthase / putrescine methyltransferase)
     - SAMS / MAT (S-adenosylmethionine synthetase)
     - ACCase-CT (acetyl-CoA carboxylase)
     - tubulins (alpha / beta tubulins)

Resolves locus IDs against TAIR10 pep FASTA (auto-fetched from Ensembl Plants if needed),
deduplicates by locus, writes essential_proteins.faa and reference_metadata.tsv,
and compiles essential_proteins.dmnd via diamond makedb.
"""

import argparse
import csv
import gzip
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from html.parser import HTMLParser
from typing import Dict, List, Optional, Set, Tuple

# Default URLs for automated remote fetching
URL_TAIR10_PEP_ENSEMBL = os.environ.get(
    "PLANTISMASH_RECRUITMENT_TAIR10_URL",
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
)
URL_SEEDGENES_HTML = os.environ.get(
    "PLANTISMASH_RECRUITMENT_SEEDGENES_URL",
    "https://seedgenes.org/GeneList.html",
)

METADATA_HEADER = ["agi_locus", "source", "family", "evidence_note"]

# Built-in curated recruited / target-pathway enzyme families
CURATED_FAMILIES_DATA: List[Tuple[str, str, str, str]] = [
    # OSC: Oxidosqualene cyclases / triterpene synthases
    ("AT1G78500", "CURATED_FAMILY", "OSC", "THAS1 (thalianol synthase 1; BGC-associated OSC)"),
    ("AT1G78510", "CURATED_FAMILY", "OSC", "MRN1 (marneral synthase 1; BGC-associated OSC)"),
    ("AT1G78950", "CURATED_FAMILY", "OSC", "PEN1 / LUP2 (oxidosqualene cyclase)"),
    ("AT1G78955", "CURATED_FAMILY", "OSC", "LUP3 (oxidosqualene cyclase)"),
    ("AT1G78960", "CURATED_FAMILY", "OSC", "LUP4 (oxidosqualene cyclase)"),
    ("AT1G78970", "CURATED_FAMILY", "OSC", "LUP1 (lupeol synthase 1)"),
    ("AT2G07050", "CURATED_FAMILY", "OSC", "CAS1 (cycloartenol synthase 1; essential primary sterol precursor)"),
    ("AT5G48010", "CURATED_FAMILY", "OSC", "BARS1 (baruol synthase 1)"),
    ("AT5G42600", "CURATED_FAMILY", "OSC", "PEN3 / LUP5 (arabidiol synthase)"),
    # CPS / KSL: Copalyl diphosphate synthase / ent-kaurene synthase (diterpene BGCs & gibberellin)
    ("AT4G02780", "CURATED_FAMILY", "CPS", "GA1 / CPS (ent-copalyl diphosphate synthase 1; gibberellin core)"),
    ("AT1G79460", "CURATED_FAMILY", "KSL", "GA2 / KS (ent-kaurene synthase 1; diterpene core)"),
    ("AT1G61680", "CURATED_FAMILY", "CPS/KSL", "TPS21 (terpene synthase 21; diterpene/sesquiterpene)"),
    ("AT1G70080", "CURATED_FAMILY", "CPS/KSL", "TPS14 (terpene synthase 14; diterpene/monoterpene)"),
    # HMGR: 3-hydroxy-3-methylglutaryl-CoA reductase (mevalonate pathway gatekeeper)
    ("AT1G76490", "CURATED_FAMILY", "HMGR", "HMG1 / HMGR1 (3-hydroxy-3-methylglutaryl-CoA reductase 1)"),
    ("AT2G17370", "CURATED_FAMILY", "HMGR", "HMG2 / HMGR2 (3-hydroxy-3-methylglutaryl-CoA reductase 2)"),
    # SPDS / PMT: Spermidine synthase / putrescine N-methyltransferase (polyamine / alkaloid recruitments)
    ("AT1G23820", "CURATED_FAMILY", "SPDS/PMT", "SPDS1 (spermidine synthase 1; primary polyamine metabolism)"),
    ("AT1G70310", "CURATED_FAMILY", "SPDS/PMT", "SPDS2 (spermidine synthase 2; primary polyamine metabolism)"),
    ("AT5G53120", "CURATED_FAMILY", "SPDS/PMT", "SPMS / SPDS3 (spermine synthase)"),
    ("AT5G19530", "CURATED_FAMILY", "SPDS/PMT", "ACL5 (thermospermine synthase ACL5)"),
    # SAMS / MAT: S-adenosylmethionine synthetase (SAM / C1 metabolism)
    ("AT1G02500", "CURATED_FAMILY", "SAMS", "SAM1 / MAT1 (S-adenosylmethionine synthetase 1)"),
    ("AT4G01850", "CURATED_FAMILY", "SAMS", "SAM2 / MAT2 (S-adenosylmethionine synthetase 2)"),
    ("AT3G17390", "CURATED_FAMILY", "SAMS", "SAM3 / MAT3 (S-adenosylmethionine synthetase 3)"),
    ("AT2G36880", "CURATED_FAMILY", "SAMS", "SAM4 / MAT4 (S-adenosylmethionine synthetase 4)"),
    # ACCase-CT: Acetyl-CoA carboxylase subunits (fatty acid / polyketide precursor generation)
    ("AT1G36160", "CURATED_FAMILY", "ACCase-CT", "ACC1 (acetyl-CoA carboxylase 1; cytosolic homomeric)"),
    ("AT1G26980", "CURATED_FAMILY", "ACCase-CT", "ACC2 (acetyl-CoA carboxylase 2; plastidic homomeric)"),
    ("AT5G16390", "CURATED_FAMILY", "ACCase-CT", "CAC1 / BCCP1 (biotin carboxyl carrier protein 1)"),
    ("AT5G35360", "CURATED_FAMILY", "ACCase-CT", "CAC2 / BCCP2 (biotin carboxyl carrier protein 2)"),
    ("AT2G38040", "CURATED_FAMILY", "ACCase-CT", "CAC3 / BC (biotin carboxylase)"),
    # Tubulin: alpha and beta tubulins (essential cytoskeletal target & resistance duplication)
    ("AT1G64740", "CURATED_FAMILY", "tubulin", "TUA1 (alpha tubulin 1)"),
    ("AT1G50010", "CURATED_FAMILY", "tubulin", "TUA2 (alpha tubulin 2)"),
    ("AT5G19770", "CURATED_FAMILY", "tubulin", "TUA3 (alpha tubulin 3)"),
    ("AT1G04820", "CURATED_FAMILY", "tubulin", "TUA4 (alpha tubulin 4)"),
    ("AT5G19780", "CURATED_FAMILY", "tubulin", "TUA5 (alpha tubulin 5)"),
    ("AT4G14960", "CURATED_FAMILY", "tubulin", "TUA6 (alpha tubulin 6)"),
    ("AT1G75780", "CURATED_FAMILY", "tubulin", "TUB1 (beta tubulin 1)"),
    ("AT5G62690", "CURATED_FAMILY", "tubulin", "TUB2 (beta tubulin 2)"),
    ("AT5G62700", "CURATED_FAMILY", "tubulin", "TUB3 (beta tubulin 3)"),
    ("AT5G44340", "CURATED_FAMILY", "tubulin", "TUB4 (beta tubulin 4)"),
    ("AT1G20010", "CURATED_FAMILY", "tubulin", "TUB5 (beta tubulin 5)"),
    ("AT5G12250", "CURATED_FAMILY", "tubulin", "TUB6 (beta tubulin 6)"),
    ("AT2G29550", "CURATED_FAMILY", "tubulin", "TUB7 (beta tubulin 7)"),
    ("AT5G23860", "CURATED_FAMILY", "tubulin", "TUB8 (beta tubulin 8)"),
    ("AT4G20890", "CURATED_FAMILY", "tubulin", "TUB9 (beta tubulin 9)"),
]


class SeedGenesHTMLParser(HTMLParser):
    """Parses SeedGenes GeneList table into structured records."""

    def __init__(self):
        super().__init__()
        self.rows: List[List[str]] = []
        self.current_row: List[str] = []
        self.in_cell = False
        self.cell_text: List[str] = []

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.current_row = []
        elif tag in ("td", "th"):
            self.in_cell = True
            self.cell_text = []

    def handle_endtag(self, tag):
        if tag in ("td", "th"):
            self.in_cell = False
            self.current_row.append(" ".join("".join(self.cell_text).split()))
        elif tag == "tr":
            if self.current_row:
                self.rows.append(self.current_row)

    def handle_data(self, data):
        if self.in_cell:
            self.cell_text.append(data)


def parse_tair_pep_text(fasta_text: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse TAIR10 protein FASTA text.
    Returns mapping: locus (uppercase, e.g. AT1G01010) -> (protein_id, sequence).
    If multiple splice variants exist, picks .1 variant or longest sequence.
    """
    seqs: Dict[str, List[Tuple[str, str]]] = {}
    current_id = None
    current_seq: List[str] = []

    for line in fasta_text.splitlines():
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
        dot_one = [v for v in variants if v[0].endswith(".1")]
        if dot_one:
            best_seqs[locus] = dot_one[0]
        else:
            best_seqs[locus] = max(variants, key=lambda v: len(v[1]))

    return best_seqs


def parse_tair_pep_fasta(fasta_path: str) -> Dict[str, Tuple[str, str]]:
    """Parse TAIR10 protein FASTA file (plain text or gzip)."""
    if fasta_path.endswith(".gz"):
        with gzip.open(fasta_path, "rt", encoding="utf-8", errors="replace") as fh:
            return parse_tair_pep_text(fh.read())
    else:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
            return parse_tair_pep_text(fh.read())


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


def try_fetch_url(url: str, timeout: int = 30) -> Optional[bytes]:
    """Attempt HTTP GET with timeout. Returns bytes or None."""
    try:
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "Mozilla/5.0 (compatible; plantismash-recruitment-builder/1.0)",
            },
        )
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return response.read()
    except Exception as e:
        logging.warning("Fetch failed for %s: %s", url, e)
        return None


def fetch_seedgenes_loci() -> List[Tuple[str, str, str, str]]:
    """Fetch and parse SeedGenes essential loci directly from SeedGenes.org."""
    logging.info("Fetching SeedGenes essential loci from %s...", URL_SEEDGENES_HTML)
    data = try_fetch_url(URL_SEEDGENES_HTML, timeout=20)
    if not data:
        logging.warning("Could not fetch SeedGenes HTML online.")
        return []

    try:
        html_text = data.decode("utf-8", errors="replace")
        parser = SeedGenesHTMLParser()
        parser.feed(html_text)

        entries = []
        for row in parser.rows:
            if not row or len(row) < 4:
                continue
            locus = row[0].strip().upper()
            if locus.startswith("AT") and len(locus) >= 9:
                symbol = row[1].strip() if len(row) > 1 else ""
                conf = row[3].strip() if len(row) > 3 else "SeedGenes"
                func = row[6].strip() if len(row) > 6 else "Essential embryo-defective gene"
                note = f"SeedGenes: {symbol} ({conf}) - {func}"
                entries.append((locus, "EMB", "EMB_essential", note))

        logging.info("Extracted %d essential loci from SeedGenes.", len(entries))
        return entries
    except Exception as e:
        logging.warning("Error parsing SeedGenes HTML: %s", e)
        return []


def run_diamond_makedb(fasta_path: str, dmnd_path: str) -> bool:
    """Run `diamond makedb --in fasta_path -d dmnd_path`."""
    diamond_bin = shutil.which("diamond")
    if not diamond_bin:
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

    db_base = dmnd_path[:-5] if dmnd_path.endswith(".dmnd") else dmnd_path
    cmd = [diamond_bin, "makedb", "--in", fasta_path, "-d", db_base]
    logging.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        logging.error("diamond makedb failed (code %d):\n%s", res.returncode, res.stderr)
        return False
    logging.info("DIAMOND database successfully created at %s", dmnd_path)
    return True


def build_recruitment_db(
    tair_pep: Optional[str] = None,
    out_dir: Optional[str] = None,
    emb_file: Optional[str] = None,
    ogee_file: Optional[str] = None,
    families_file: Optional[str] = None,
    fetch_missing: bool = True,
) -> bool:
    """Main automated build workflow."""
    if not out_dir:
        out_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "antismash", "generic_modules", "recruitment_miner", "data")
        )
    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse TAIR10 pep (use local file or auto-download from Ensembl Plants)
    tair_seqs: Dict[str, Tuple[str, str]] = {}
    if tair_pep and os.path.isfile(tair_pep):
        logging.info("Parsing TAIR10 peptide sequences from %s...", tair_pep)
        tair_seqs = parse_tair_pep_fasta(tair_pep)
    elif fetch_missing:
        logging.info("TAIR10 pep not provided locally; fetching Arabidopsis TAIR10 proteome from Ensembl Plants: %s...", URL_TAIR10_PEP_ENSEMBL)
        data = try_fetch_url(URL_TAIR10_PEP_ENSEMBL, timeout=120)
        if data:
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
                fasta_text = gz.read().decode("utf-8", errors="replace")
            tair_seqs = parse_tair_pep_text(fasta_text)
        else:
            logging.error("Could not download TAIR10 pep from Ensembl Plants.")
            return False
    else:
        logging.error("TAIR10 pep FASTA file not found and remote fetching disabled: %s", tair_pep)
        return False

    logging.info("Parsed %d unique AGI loci from TAIR10 reference proteome.", len(tair_seqs))

    # 2. Collect candidate entries
    entries: List[Tuple[str, str, str, str]] = []

    # Built-in curated enzyme families
    entries.extend(CURATED_FAMILIES_DATA)

    if families_file and os.path.isfile(families_file):
        logging.info("Parsing extra curated families from %s...", families_file)
        entries.extend(parse_families_file(families_file))

    if emb_file and os.path.isfile(emb_file):
        logging.info("Parsing EMB loci from %s...", emb_file)
        entries.extend(parse_emb_file(emb_file))
    elif fetch_missing:
        # Auto-fetch SeedGenes loci online
        seedgenes_entries = fetch_seedgenes_loci()
        entries.extend(seedgenes_entries)

    if ogee_file and os.path.isfile(ogee_file):
        logging.info("Parsing OGEE entries from %s...", ogee_file)
        entries.extend(parse_ogee_file(ogee_file))

    if not entries:
        logging.warning("No reference entries collected! Check input files or network.")
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
                deduped[locus] = (source, family, f"{curr_note}; {note}" if curr_note else note)
            elif source == curr_src:
                if family and family != curr_fam and curr_fam == "EMB_essential":
                    deduped[locus] = (source, family, f"{curr_note}; {note}" if curr_note else note)

    logging.info("Total unique reference loci to process: %d", len(deduped))

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
        description="Automated reference DB builder for recruitment_miner (EMB/SeedGenes, OGEE, and Curated Families)."
    )
    parser.add_argument("--tair-pep", default=None, help="TAIR10 pep FASTA path (optional; auto-downloads if omitted)")
    parser.add_argument("--emb-loci", dest="emb_file", default=None, help="EMB / SeedGenes locus file (optional)")
    parser.add_argument("--ogee", dest="ogee_file", default=None, help="OGEE essential gene TSV (optional)")
    parser.add_argument("--families", dest="families_file", default=None, help="Curated recruited families TSV (optional)")
    parser.add_argument("--out-dir", default=default_out, help="Output directory for essential DB")
    parser.add_argument("--no-fetch", dest="fetch", action="store_false", default=True, help="Disable online fetching")

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
        sys.exit(1)
    else:
        logging.info("Recruitment DB build completed successfully.")


if __name__ == "__main__":
    main()

