# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
from substitution_matrices import *
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    score, aligned_seq1, aligned_seq2 = nw.align(seq1, seq2)
    #ADD MATRIX ASSERT STATEMENTS
    assert len(aligned_seq1) == len(aligned_seq2)
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    score, aligned_seq3, aligned_seq4 = nw.align(seq3, seq4)
    assert(aligned_seq3 == "MAVHQLIRRP")
    assert(aligned_seq4 == "M---QLIRHP")
    assert(score == 17)
    assert len(aligned_seq3) == len(aligned_seq4)




