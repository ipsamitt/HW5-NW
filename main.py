# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    scores = []
    alignments = []
    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    gg_score, gg_align, hs_gg_align = nw.align(gg_seq, hs_seq)
    scores.append(gg_score)
    alignments.append("Gallus gallus")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    mm_score, mm_align, hs_mm_align = nw.align(mm_seq, hs_seq)
    scores.append(mm_score)
    alignments.append("Mus musculus")


    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    br_score, br_align, hs_br_align = nw.align(br_seq, hs_seq)  
    scores.append(br_score)
    alignments.append("Balaeniceps rex")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    tt_score, tt_align, hs_tt_align = nw.align(tt_seq, hs_seq)
    scores.append(tt_score)
    alignments.append("Tursiops_truncatus")

    rank_scores = (np.argsort(scores))[::-1] 

    for i in rank_scores:
        print(alignments[i])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for i in rank_scores:
        print(alignments[i] + " " + str(scores[i]))
    

if __name__ == "__main__":
    main()
