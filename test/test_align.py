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
    #calculated by hand
    true_gapA = np.array([[-10, -float('inf'), -float('inf'), -float('inf')], 
                          [-11, -float('inf'), -float('inf'), -float('inf')], 
                          [-12, -6, -22, -24], 
                          [-13, -7, -7, -19], 
                          [-14, -8, -8, -6]])
    assert np.allclose(nw._gapA_matrix, true_gapA, atol=1e-6)
    true_gapB = np.array([[-10, -11, -12, -13],
                           [-float('inf'), -float('inf'), -6, -7], 
                           [-float('inf'), -float('inf'), -23, -7],
                             [-float('inf'), -float('inf'), -23, -12], 
                             [-float('inf'), -float('inf'), -25, -17]])
    assert np.allclose(nw._gapB_matrix, true_gapB, atol=1e-6)
    true_align = np.array([[0, -float('inf'), -float('inf'),-float('inf')],
                            [-float('inf'), 5, -11, -13], 
                            [-float('inf'), -12, 4, -8], 
                            [-float('inf'), -12, -1, 5], 
                            [-float('inf'), -14, -6, 4]])
    assert np.allclose(nw._align_matrix, true_align, atol=1e-6)

def test_nw_backtrace():
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    #check alignments and scores
    score, aligned_seq3, aligned_seq4 = nw.align(seq3, seq4)
    assert(aligned_seq3 == "MAVHQLIRRP")
    assert(aligned_seq4 == "MQ---LIRHP")
    assert(score == 17)



