#!/usr/bin/env python3
"""Build plants_extended model dir: extract curated dbCAN + Pfam HMM profiles,
rewrite NAME/DESC headers to plantismash conventions (prefixed ids), append
hmmdetails.txt entries."""

import gzip
import sys
from pathlib import Path

EXT = Path(sys.argv[1])  # plants_extended dir

# family -> new profile id, description
DBCAN_SELECT = {
    "GT1.hmm": ("dbcan_GT1", "UDP-glycosyltransferase family (CAZy GT1)"),
    "GT2_Cellulose_synt.hmm": (
        "dbcan_GT2_cellulose_synt",
        "Cellulose synthase-like glycosyltransferase (CAZy GT2)",
    ),
    "GH1.hmm": ("dbcan_GH1", "Beta-glucosidase (CAZy GH1)"),
    "GH3.hmm": (
        "dbcan_GH3",
        "Glycoside hydrolase family 3 (beta-glucosidase/xylosidase)",
    ),
    "AA3.hmm": ("dbcan_AA3", "GMC-family oxidoreductase (CAZy AA3)"),
    "AA7.hmm": ("dbcan_AA7", "Glucooligosaccharide oxidase, FAD-dependent (CAZy AA7)"),
    "CE1.hmm": ("dbcan_CE1", "Acetylxylan esterase family (CAZy CE1)"),
}
PFAM_SELECT = {
    "PF00685": ("pfam_Sulfotransferase", "Sulfotransferase domain (Pfam PF00685)"),
    "PF00043": (
        "pfam_GST_C",
        "Glutathione S-transferase C-terminal domain (Pfam PF00043)",
    ),
    "PF02798": (
        "pfam_GST_N",
        "Glutathione S-transferase N-terminal domain (Pfam PF02798)",
    ),
    # Primary-metabolism suppressors: compete with CAZyme labels via overlap groups;
    # when these win on a gene, the gene is excluded from BGC rule counting.
    "PF03552": (
        "pfam_CESA",
        "Cellulose synthase, primary metabolism suppressor (Pfam PF03552)",
    ),
    "PF00534": (
        "pfam_SUSY",
        "Glycosyl transferases group 1 incl. sucrose synthase, primary metabolism suppressor (Pfam PF00534)",
    ),
    "PF01095": (
        "pfam_PME",
        "Pectin methylesterase, primary metabolism suppressor (Pfam PF01095)",
    ),
}

HEADER_KEYS_SKIP = {"NAME", "ACC", "DESC"}  # lines we rewrite


def parse_blocks(path, opener=open):
    """Yield HMM text blocks terminated by '//'."""
    buf = []
    with opener(path) as handle:
        for line in handle:
            buf.append(line)
            if line.rstrip() == "//":
                yield "".join(buf)
                buf = []


def block_name(block):
    for line in block.splitlines():
        if line.startswith("NAME"):
            return line.split(None, 1)[1].strip()
    return None


def rewrite_block(block, new_name, desc):
    out = []
    for line in block.splitlines(keepends=True):
        if line.startswith("NAME"):
            out.append(f"NAME  {new_name}\n")
        elif line.startswith("ACC"):
            out.append(line)
        elif line.startswith("DESC"):
            out.append(f"DESC  {desc}\n")
        else:
            out.append(line)
    return "".join(out)


def extract(src, wanted_by_name, prefix_rewriter):
    found = {}
    for block in parse_blocks(src):
        name = block_name(block)
        if name in wanted_by_name:
            new_id, desc = wanted_by_name[name]
            target = EXT / f"{new_id}.hmm"
            target.write_text(prefix_rewriter(new_id, desc, block))
            found[name] = new_id
            print(f"  extracted {name} -> {target.name}")
    missing = set(wanted_by_name) - set(found)
    if missing:
        raise SystemExit(f"ERROR: models not found in source: {missing}")
    return found


def main():
    details = []

    # --- dbCAN ---
    print("== dbCAN ==")

    def rw_dbcan(new_id, desc, block):
        # strip .hmm from NAME handled by rewrite_block; keep ACC as original family
        return rewrite_block(block, new_id, desc)

    extract(
        "/tmp/bgc_meta/dbCAN.hmm", {k: v for k, v in DBCAN_SELECT.items()}, rw_dbcan
    )

    # --- Pfam ---
    print("== Pfam ==")

    def rw_pfam(new_id, desc, block):
        return rewrite_block(block, new_id, desc)

    # stream through gzip directly, extracting blocks by accession
    pfam_by_acc = {acc: v for acc, v in PFAM_SELECT.items()}
    found_pfam = {}
    buf = []
    with gzip.open("/tmp/bgc_meta/Pfam-A.hmm.gz", "rt") as fh:
        for line in fh:
            buf.append(line)
            if line.rstrip() == "//":
                block = "".join(buf)
                buf = []
                acc = None
                for l in block.splitlines():
                    if l.startswith("ACC"):
                        acc = l.split(None, 1)[1].strip().split(".")[0]
                        break
                if acc in pfam_by_acc:
                    new_id, desc = pfam_by_acc[acc]
                    (EXT / f"{new_id}.hmm").write_text(rw_pfam(new_id, desc, block))
                    found_pfam[acc] = new_id
                    print(f"  extracted {acc} -> {new_id}.hmm")
    missing = set(pfam_by_acc) - set(found_pfam)
    if missing:
        raise SystemExit(f"ERROR: Pfam accessions not found: {missing}")

    # --- hmmdetails.txt append (idempotent) ---
    details_path = EXT / "hmmdetails.txt"
    existing = {
        line.split("\t")[0]
        for line in details_path.read_text().splitlines()
        if line.strip()
    }
    new_lines = [
        f"{pid}\t{desc}\t-1\t{pid}.hmm"
        for pid, desc in [v for v in DBCAN_SELECT.values()]
        + [v for v in PFAM_SELECT.values()]
        if pid not in existing
    ]
    with details_path.open("a") as fh:
        for line in new_lines:
            fh.write(line + "\n")
    print(f"appended {len(new_lines)} hmmdetails entries")

    # --- cluster_rules.txt extension ---
    # Reuse the baseline general-pool bracket verbatim from the 'plant' rule,
    # then append CAZyme-anchored rules. Idempotent on rule name.
    rules_path = EXT / "cluster_rules.txt"
    rules_text = rules_path.read_text()
    plant_line = next(l for l in rules_text.splitlines() if l.startswith("plant\t"))
    pool = plant_line.partition("[")[2].partition("]")[0]  # general SM enzyme list

    glyco_pool = pool + ",dbcan_GT1,dbcan_GT2_cellulose_synt,dbcan_GH1,dbcan_GH3"
    cazyme_pool = glyco_pool + ",dbcan_AA3,dbcan_AA7,dbcan_CE1"

    new_rules = []
    if "\nglycoside\t" not in rules_text:
        # glycosylated-metabolite pattern: >=4 pool genes nearby, AND a
        # glycosyltransferase-signature gene AND a beta-glucosidase-signature gene
        sig = "dbcan_GT1/dbcan_GT2_cellulose_synt/UDPGT/UDPGT_2/Glycos_transf_1/Glycos_transf_2"
        ess = "dbcan_GH1/dbcan_GH3"
        new_rules.append(
            f"glycoside\tminimum(4,[{glyco_pool}],[{sig},{ess}])\t5\t1"
        )
    if "\ncazyme_rich\t" not in rules_text:
        # CAZyme-enriched oxidative/conjugation tailoring regions (experimental)
        sig = "dbcan_AA7/dbcan_AA3/dbcan_CE1/pfam_Sulfotransferase/pfam_GST_C"
        new_rules.append(
            f"cazyme_rich\tminimum(4,[{cazyme_pool}],[{sig}])\t5\t1"
        )
    if new_rules:
        with rules_path.open("a") as fh:
            for line in new_rules:
                fh.write(line + "\n")
        print(f"appended {len(new_rules)} cluster rules")
    else:
        print("cluster rules already present")

if __name__ == "__main__":
    main()
