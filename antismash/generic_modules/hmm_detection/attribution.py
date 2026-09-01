import os
import logging
_hits = []
_overlaps = []
_clusters = []
_enabled = True
_output_dir = None
def init(output_dir, enabled=True):
    global _hits, _overlaps, _clusters, _enabled, _output_dir
    _hits = []
    _overlaps = []
    _clusters = []
    _enabled = enabled
    _output_dir = output_dir
def is_enabled():
    return _enabled
def record_hit(contig, gene_id, profile, bits, evalue, aln_len, model_len, coverage, decision, reason):
    if not _enabled:
        return
    _hits.append([contig, gene_id, profile, bits, evalue, aln_len, model_len, coverage, decision, reason])
def record_overlap(contig, gene_id, group_members, winner, winner_bits, loser, loser_bits, reason):
    if not _enabled:
        return
    _overlaps.append([contig, gene_id, group_members, winner, winner_bits, loser, loser_bits, reason])
def record_cluster(cluster_num, contig, start, end, type_, rules, n_genes, n_cores, marginal_quorum, quorum_detail, label_losses):
    if not _enabled:
        return
    _clusters.append([cluster_num, contig, start, end, type_, rules, n_genes, n_cores, marginal_quorum, quorum_detail, label_losses])
def write():
    if not _enabled or _output_dir is None:
        return
    outdir = os.path.join(_output_dir, "attribution")
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError:
        pass
    try:
        with open(os.path.join(outdir, "hits.tsv"), "w") as fh:
            fh.write("contig\tgene_id\tprofile\tbits\tevalue\taln_len\tmodel_len\tcoverage\tdecision\treason\n")
            for r in _hits:
                fh.write("\t".join(map(str, r)) + "\n")
    except OSError as e:
        logging.warning("attribution hits write failed: %s", e)
    try:
        with open(os.path.join(outdir, "overlap_decisions.tsv"), "w") as fh:
            fh.write("contig\tgene_id\tgroup_members\twinner\twinner_bits\tloser\tloser_bits\treason\n")
            for r in _overlaps:
                fh.write("\t".join(map(str, r)) + "\n")
    except OSError as e:
        logging.warning("attribution overlaps write failed: %s", e)
    try:
        with open(os.path.join(outdir, "cluster_summary.tsv"), "w") as fh:
            fh.write("cluster_num\tcontig\tstart\tend\ttype\trules\tn_genes\tn_cores\tmarginal_quorum\tquorum_detail\tlabel_losses\n")
            for r in _clusters:
                fh.write("\t".join(map(str, r)) + "\n")
    except OSError as e:
        logging.warning("attribution clusters write failed: %s", e)
def clear():
    global _hits, _overlaps, _clusters
    _hits = []
    _overlaps = []
    _clusters = []
