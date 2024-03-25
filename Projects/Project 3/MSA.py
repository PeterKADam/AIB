import alignment
from itertools import combinations
import numpy as np
import sys

score_matrix = {
    "A": {"A": 0, "C": 5, "G": 2, "T": 5},
    "C": {"A": 5, "C": 0, "G": 5, "T": 2},
    "G": {"A": 2, "C": 5, "G": 0, "T": 5},
    "T": {"A": 5, "C": 2, "G": 5, "T": 0},
}
gap_penalty = 5

short_seq = {
    "seq1": "GTTCCGAAAGGCTAGCGCTAGGCGCC",
    "seq2": "ATGGATTTATCTGCTCTTCG",
    "seq3": "TGCATGCTGAAACTTCTCAACCA",
}
short_score = 198
long_seq = {
    "seq1": "GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGATAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTCTCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTG",
    "seq2": "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAACGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA",
    "seq3": "CGCTGGTGCAACTCGAAGACCTATCTCCTTCCCGGGGGGGCTTCTCCGGCATTTAGGCCTCGGCGTTTGGAAGTACGGAGGTTTTTCTCGGAAGAAAGTTCACTGGAAGTGGAAGAAATGGATTTATCTGCTGTTCGAATTCAAGAAGTACAAAATGTCCTTCATGCTATGCAGAAAATCTTGGAGTGTCCAATCTGTTT",
}


# exact MSA


def traceback_arrows(i, j, k, seq1_base, seq2_base, seq3_base, D):

    node_score = D[i][j][k]

    up = D[i - 1][j][k]
    left = D[i][j - 1][k]
    depth = D[i][j][k - 1]

    diag_up = D[i - 1][j - 1][k]
    diag_left = D[i - 1][j][k - 1]
    diag_depth = D[i][j - 1][k - 1]

    diag_all = D[i - 1][j - 1][k - 1]

    match_score_all = (
        score_matrix[seq1_base][seq2_base] + score_matrix[seq1_base][seq3_base] + score_matrix[seq2_base][seq3_base]
    )

    match_score_up = score_matrix[seq1_base][seq2_base]
    match_score_left = score_matrix[seq1_base][seq3_base]
    match_score_depth = score_matrix[seq2_base][seq3_base]

    if node_score == diag_all + match_score_all:
        return "diag_all"
    elif node_score == diag_up + match_score_up:
        return "diag_up"
    elif node_score == diag_left + match_score_left:
        return "diag_left"
    elif node_score == diag_depth + match_score_depth:
        return "diag_depth"
    elif node_score == up + (2 * gap_penalty):
        return "up"
    elif node_score == left + (2 * gap_penalty):
        return "left"
    elif node_score == depth + (2 * gap_penalty):
        return "depth"

    sys.exit("Error in traceback_arrows: no match found.")


def traceback(D, seq_dict):
    seq1, seq2, seq3 = seq_dict["seq1"], seq_dict["seq2"], seq_dict["seq3"]
    i, j, k = len(seq1), len(seq2), len(seq3)
    seq1_align, seq2_align, seq3_align = "", "", ""
    while i > 0 or j > 0 or k > 0:
        arrows = traceback_arrows(i, j, k, seq1[i - 1], seq2[j - 1], seq3[k - 1], D)
        if arrows == "diag_all":
            seq1_align += seq1[i - 1]
            seq2_align += seq2[j - 1]
            seq3_align += seq3[k - 1]
            i -= 1
            j -= 1
            k -= 1
        elif arrows == "diag_up":
            seq1_align += seq1[i - 1]
            seq2_align += seq2[j - 1]
            seq3_align += "-"
            i -= 1
            j -= 1
        elif arrows == "diag_left":
            seq1_align += seq1[i - 1]
            seq2_align += "-"
            seq3_align += seq3[k - 1]
            i -= 1
            k -= 1
        elif arrows == "diag_depth":
            seq1_align += "-"
            seq2_align += seq2[j - 1]
            seq3_align += seq3[k - 1]
            j -= 1
            k -= 1
        elif arrows == "up":
            seq1_align += seq1[i - 1]
            seq2_align += seq2[j - 1]
            seq3_align += "-"
            i -= 1
            j -= 1
        elif arrows == "left":
            seq1_align += seq1[i - 1]
            seq2_align += "-"
            seq3_align += seq3[k - 1]
            i -= 1
            k -= 1
        elif arrows == "depth":
            seq1_align += "-"
            seq2_align += seq2[j - 1]
            seq3_align += seq3[k - 1]
            j -= 1
            k -= 1
        # print(f"{arrows}")
        # print(f"{seq1_align}, {seq2_align}, {seq3_align}")
    return seq1_align[::-1], seq2_align[::-1], seq3_align[::-1]


def D_calc(seq_dict):
    seq1, seq2, seq3 = seq_dict["seq1"], seq_dict["seq2"], seq_dict["seq3"]

    # D = empty_matrix(seq_dict)
    D = np.full((len(seq1) + 1, len(seq2) + 1, len(seq3) + 1), None)
    # print(f"D:{len(D)}{len(D[0])}{len(D[[0][0]])}")

    # it really helps to cache the pairwise alignments, because they are used a lot, theyre probably gonna be cpu cached anyway, but whatever, now theyre explicit.
    # and [i],[j] is the same as [:i], [:j] in the alignment matrix (eg, you dont have to calculate the smaller alignment for every step.)
    # these have to have reversed indexes, because reasons
    alignmentscache_ij = alignment.LinearGlobalAlignment(seq2, seq1, score_matrix, gap_penalty).matrix
    alignmentscache_ik = alignment.LinearGlobalAlignment(seq3, seq1, score_matrix, gap_penalty).matrix
    alignmentscache_jk = alignment.LinearGlobalAlignment(seq3, seq2, score_matrix, gap_penalty).matrix

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            D[i][j][0] = alignmentscache_ij[i][j] + (i + j) * gap_penalty
            # print(D[i][j][0])
    for i in range(len(seq1)):
        for k in range(len(seq3)):
            D[i][0][k] = alignmentscache_ik[i][k] + (i + k) * gap_penalty
    for j in range(len(seq2)):
        for k in range(len(seq3)):
            D[0][j][k] = alignmentscache_jk[j][k] + (k + j) * gap_penalty

    # print(f"D:{D[0][0]}")

    for i in range(0, len(seq1) + 1):
        for j in range(0, len(seq2) + 1):
            for k in range(0, len(seq3) + 1):
                # basei = seq1[i]
                # basej = seq2[j]
                # basek = seq3[k]

                # This is technically reverse indexing ie [-1] on the first pass, but the scores arent used at any [0] so yolo
                scores = [
                    score_matrix[seq1[i - 1]][seq2[j - 1]],
                    score_matrix[seq1[i - 1]][seq3[k - 1]],
                    score_matrix[seq2[j - 1]][seq3[k - 1]],
                ]

                v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")

                if i == 0 and j == 0 and k == 0:
                    v0 = 0

                if i > 0 and j > 0 and k > 0:
                    v1 = D[i - 1][j - 1][k - 1] + scores[0] + scores[1] + scores[2]
                if i > 0 and j > 0 and k >= 0:
                    v2 = D[i - 1][j - 1][k] + scores[0] + (2 * gap_penalty)
                if i > 0 and j >= 0 and k > 0:
                    v3 = D[i - 1][j][k - 1] + scores[1] + (2 * gap_penalty)
                if i >= 0 and j > 0 and k > 0:
                    v4 = D[i][j - 1][k - 1] + scores[2] + (2 * gap_penalty)
                if i > 0 and j >= 0 and k >= 0:
                    v5 = D[i - 1][j][k] + (2 * gap_penalty)
                if i >= 0 and j > 0 and k >= 0:
                    v6 = D[i][j - 1][k] + (2 * gap_penalty)
                if i >= 0 and j >= 0 and k > 0:
                    v7 = D[i][j][k - 1] + (2 * gap_penalty)
                D[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)
                # print(f"D:{D[i][j][k]}\ni:{i}\nj:{j}\nk:{k}\n")

        if (i / len(seq1)) * 100 % 5 == 0 and i != 0:
            print(f"{(i / len(seq1)) * 100}%")

    return D


##approximate MSA


def find_center_string(seqs, score, gap):
    # find the center string by comparing all the sequences to each other
    min_score = float("inf")
    center_key = None
    for key, seq in seqs.items():
        score_sum = 0
        for key2, seq2 in seqs.items():
            if key != key2:
                score_sum += alignment.LinearGlobalAlignment(seq, seq2, score, gap).get_max()
        if score_sum < min_score:
            min_score = score_sum
            center_key = key
    return center_key


def extend_approx_MSA(M, Alignment):

    return M


def MSA(seqs, score, gap):
    center_key = find_center_string(seqs, score, gap)
    center_seq = seqs.pop(center_key)
    alignments = []
    for key, seq in seqs.items():
        all_alignments = alignments.append(
            alignment.LinearGlobalAlignment(center_seq, seq, score, gap).get_alignments()[0]
        )
    M = []




# print(find_center_string(long_seq, score_matrix, gap_penalty))
MSA(short_seq, score_matrix, gap_penalty)
