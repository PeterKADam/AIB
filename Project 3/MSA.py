import alignment
from itertools import combinations
import numpy as np

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


def get_Pairwisealignments(seqs, gap_penalty, score_matrix):
    alignments = {}
    for seq1, seq2 in combinations(seqs.keys(), 2):
        alignments[(seq1, seq2)] = alignment.LinearGlobalAlignment(
            seqs[seq1], seqs[seq2], score_matrix, gap_penalty
        ).alignments[0][0]
    return alignments


def empty_matrix(seqs):
    return [
        [
            [None for _ in range(len(seqs["seq1"]) + 1)]
            for _ in range(len(seqs["seq2"]) + 1)
        ]
        for _ in range(len(seqs["seq3"]) + 1)
    ]


def D_calc(seq_dict):
    seq1, seq2, seq3 = seq_dict["seq1"], seq_dict["seq2"], seq_dict["seq3"]

    # D = empty_matrix(seq_dict)
    D = np.full((len(seq1), len(seq2), len(seq3)), None)
    # print(f"D:{len(D)}{len(D[0])}{len(D[[0][0]])}")
    len_i = len(seq1)
    len_j = len(seq2)
    len_k = len(seq3)

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            for k in range(len(seq3)):
                scores = [
                    score_matrix[seq1[i]][seq2[j]],
                    score_matrix[seq1[i]][seq3[k]],
                    score_matrix[seq2[j]][seq3[k]],
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

    return D[-1][-1][-1]


print(D_calc(short_seq))
# print(D_calc(long_seq))
