#!/usr/bin/env python3
import re
import sys
sys.path.insert(0, ".")
DOM_RE = re.compile(r'([A-Za-z0-9_\.\-]+) \(E-value')
def _mk_hit(qid, bits, hstart, hend):
    class H: pass
    h=H(); h.query_id=qid; h.bitscore=bits; h.hit_start=hstart; h.hit_end=hend; h.hit_id="g1"; h.evalue=1e-10; h._aln_len=100; h._model_len=300; h._coverage=0.33; return h
def test_length_gate():
    cov=0.25
    for ml, aln, expect_keep in [(93, 30, True), (76, 20, True), (722, 50, False), (722, 200, True), (298, 74, True)]:
        cov_ok = (ml < 300) or (aln/ml >= cov)
        assert cov_ok == expect_keep, "gate %s %s got %s"%(ml,aln,cov_ok)
    print("length_gate ok")
def test_floor_override():
    def keep(bits, cutoff, floor):
        return (bits > cutoff) and (floor is None or bits >= floor)
    assert keep(30, 25, 26)==True
    assert keep(25.5, 25, 26)==False  # AND: floor 26 tightens cutoff 25
    assert keep(60, 25, 60)==True
    assert keep(20, 25, None)==False
    assert keep(26, 25, 27)==False  # AND: 26>25 but 26<27 dropped
    assert keep(52.1, 25, 52)==True  # GH3 edge from calibration
    print("floor_override ok")
def test_suppressor_margin():
    def suppressor_wins(winner_bits, runner_bits, margin=0.15):
        return winner_bits >= (1+margin)*runner_bits
    assert suppressor_wins(110, 96.5, 0.15) == False
    assert suppressor_wins(120, 96.5, 0.15) == True
    assert suppressor_wins(111, 96.5, 0.15) == True
    assert suppressor_wins(110.9, 96.5, 0.15) == False
    print("suppressor_margin ok")
def test_marginal_quorum():
    def is_marginal(n_genes, min_number, bits_list):
        return n_genes==min_number and len(bits_list)==min_number and all(b<50 for b in bits_list)
    assert is_marginal(4,4,[30,40,45,49])==True
    assert is_marginal(4,4,[30,40,45,55])==False
    assert is_marginal(5,4,[30,40,45,49])==False
    assert is_marginal(4,4,[])==False
    print("marginal_quorum ok")
def test_sec_met_regex():
    s="pfam_GST_C (E-value: 1e-10, bitscore: 42.5, seeds: 10, aln/model=80/93, cov=0.86)"
    assert DOM_RE.search(s) is not None
    assert DOM_RE.search(s).group(1)=="pfam_GST_C"
    s2="dbcan_AA7 (E-value: 1e-20, bitscore: 100, seeds: 5, aln/model=200/458, cov=0.44)"
    assert DOM_RE.search(s2) is not None
    print("sec_met_regex ok")
if __name__=="__main__":
    test_length_gate(); test_floor_override(); test_suppressor_margin(); test_marginal_quorum(); test_sec_met_regex()
    print("all unit tests passed")
