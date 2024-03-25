import alignment
from itertools import combinations
import numpy as np
import sys
from Bio import SeqIO
from Bio.Seq import Seq

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
    elif node_score == diag_up + match_score_up + (2 * gap_penalty):
        return "diag_up"
    elif node_score == diag_left + match_score_left + (2 * gap_penalty):
        return "diag_left"
    elif node_score == diag_depth + match_score_depth + (2 * gap_penalty):
        return "diag_depth"
    elif node_score == up + (2 * gap_penalty):
        return "up"
    elif node_score == left + (2 * gap_penalty):
        return "left"
    elif node_score == depth + (2 * gap_penalty):
        return "depth"

    sys.exit("Error in traceback_arrows: no match found.")


def traceback(D, seq_dict):
    seq1, seq2, seq3 = seq_dict.values()
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


def D_calc(seq_dict, score_matrix, gap_penalty):
    seq1, seq2, seq3 = seq_dict.values()

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
            print("\033[1A", end="\x1b[2K")
            print(f"{(i / len(seq1)) * 100}%")

    return D


def sp_exact_3(seq_dict, score, gap, alignment: bool = False):
    D = D_calc(seq_dict)
    if not alignment:
        return D[-1][-1][-1]
    else:
        return traceback(D_calc(seq_dict), seq_dict)


##approximate MSA


def find_center_key(seqs, score, gap):
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


def MSA2(seqs, score, gap):
    center_key = find_center_key(seqs, score, gap)
    center_seq = seqs.pop(center_key)
    alignments = []  # a list is exceptionally terrible to use here.

    for key, seq in seqs.items():
        all_alignments = alignments.append(
            alignment.LinearGlobalAlignment(center_seq, seq, score, gap).get_alignments()[0]
        )

    M = []


def strings_to_lists(alignment):
    # Use zip function to iterate over corresponding characters of each string
    # Then use list comprehension to create the list of lists
    return [list(chars) for chars in zip(*alignment)]


def lists_to_strings(list_of_lists):
    # Use zip function to iterate over corresponding elements of each list
    # Then join the characters together and append to a new list
    return ["".join(chars) for chars in zip(*list_of_lists)]


def MSA(seqs, score, gap):
    center_key = find_center_key(seqs, score, gap)
    center_seq = seqs.pop(center_key)
    alignments = []  # a list is exceptionally terrible to use here.

    for key, seq in seqs.items():
        if center_seq != seq:
            A1 = alignment.LinearGlobalAlignment(seq, center_seq, score, gap).get_alignments()[0]
            A2 = strings_to_lists(A1)
            alignments.append(A2)

    M = alignments[0]
    MA = []

    for a in range(1, len(alignments)):
        i = 0
        j = 0
        A = alignments[a]

        while i < len(M) and j < len(A):

            # Invariant: (1) MA is a valid merge of all columns before column i in M
            # and all columns before column in A, and (2) the first row of M and A up
            # to (but not including) column i and j respectively is the same string
            # if gaps are removed.

            if M[i][0] == "-" and A[j][0] == "-":
                # Case 1: The next column in MA is column i in M extended with the second symbol
                # in column j in A.
                M[i].append(A[j][1])
                MA.append(M[i])
                i = i + 1
                j = j + 1

            elif M[i][0] == "-" and A[j][0] != "-":
                # Case 2: A[j][0] is a character, so the second symbol in column j in A, A[j][1],
                # must be in the column of MA that is the column in M where the first symbol corresponds
                # to A[j][0]. By the invariant, this column in M is the next column in M, where the first
                # symbol is a character, so we just moved forward in M until we find this column.
                M[i].append("-")
                MA.append(M[i])
                i = i + 1

            elif M[i][0] != "-" and A[j][0] == "-":
                # Case 3: M[i][0] is a character, so column i in M must be in the column of MA that also
                # contains the second symbol from the column in A, where the first symbol is the character
                # corresponding to M[i][0]. By the invariant, this column in A is the next column in A,
                # where the first symbol is a character, so we just add columns from A to MA until we
                # find this column.
                c = ["-"] * len(M[i])
                c.append(A[j][1])
                MA.append(c)
                j = j + 1

            elif M[i][0] != "-" and A[j][0] != "-":
                # Case 4: By the invariant the characters M[i][0] and A[j][0] are at the same position
                # in the string spelled by the row of M and A if gaps are removed. The next column in
                # MA is thus column i in M extended with the second symbol in column j in A.
                M[i].append(A[j][1])
                MA.append(M[i])
                i = i + 1
                j = j + 1

        if i < len(M):
            # add the remaining coloumns of M to MA
            while i < len(M):
                MA.append(M[i].append("-"))
                i = i + 1

        if j < len(A):
            # add the remaining columns of A to MA
            k = len(MA[-1])
            while j < len(A):
                c = ["-"] * (k - 1)
                c.append(A[j][1])
                MA.append(c)
                j = j + 1

    # MA2 = lists_to_strings(MA)
    return MA


# print(find_center_string(long_seq, score_matrix, gap_penalty))
MSA(short_seq, score_matrix, gap_penalty)

# print(traceback(D_calc(short_seq), short_seq))

### load files for tests


def load_match_scores(path):
    MATCH_SCORES = {}

    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                values = line.split()
                key = values.pop(0)
                MATCH_SCORES[key] = {k: int(v) for k, v in zip(["A", "C", "G", "T"], values)}
    return MATCH_SCORES


def load_sequences(filepath):
    sequences = {}
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


### run tests

### experiment 1


def experiment_1():
    ex1_seqs = load_sequences(r"Projects\Project 3\data\brca1-testseqs.fasta")
    ex1_seqs = dict(list(ex1_seqs.items())[:3])

    ex1_optimal = sp_exact_3(ex1_seqs, score_matrix, gap_penalty, alignment=False)

    ex1_alignment = sp_exact_3(ex1_seqs, score_matrix, gap_penalty, alignment=True)

    print(f"Experiment 1 optimal score: {ex1_optimal}")
    for i, alignment in enumerate(ex1_alignment, start=1):
        print(f"{alignment}")


# experiment_1()


## experiment 2
def experiment2_():
    ex2_seqs = load_sequences(r"Projects\Project 3\data\brca1-testseqs.fasta")
    ex2_seqs = dict(list(ex2_seqs.items())[:5])

    ex2_center_key = find_center_key(ex2_seqs, score_matrix, gap_penalty)
    ex2_optimal = MSA(ex2_seqs, score_matrix, gap_penalty)

    print(f"Experiment 2 optimal score: {ex2_optimal}\n Center sequence: {ex2_center_key}")


# experiment2_()


### experiment 3
def experiment3():
    import matplotlib.pyplot as plt

    # Experiment 3
    ratios = []
    lengths = []

    for i in range(10, 210, 10):
        ex3_seqs = load_sequences(f"Projects/Project 3/data/testseqs_{i}_3.fasta")
        ex3_seqs = dict(list(ex3_seqs.items())[:3])

        ex3_optimal_exact = sp_exact_3(ex3_seqs, score_matrix, gap_penalty)
        ex3_optimal_approx = sp_approx(ex3_seqs, score_matrix, gap_penalty)

        ratio = ex3_optimal_approx / ex3_optimal_exact
        ratios.append(ratio)
        lengths.append(i)

        print(f"Experiment 3 sequence length {i} ratio: {ratio}")

    # Plot the ratios
    plt.plot(lengths, ratios)
    plt.xlabel("Sequence Length")
    plt.ylabel("Approximation Ratio")
    plt.title("Approximation Ratio of sp_approx to sp_exact_3")
    plt.show()


def presentation():
    seq_sets = [
        ["brca1_bos_taurus", "brca1_canis_lupus", "brca1_gallus_gallus"],
        ["brca1_bos_taurus", "brca1_canis_lupus", "brca1_gallus_gallus", "brca1_homo_sapiens"],
        ["brca1_bos_taurus", "brca1_canis_lupus", "brca1_gallus_gallus", "brca1_homo_sapiens", "brca1_macaca_mulatta"],
        [
            "brca1_bos_taurus",
            "brca1_canis_lupus",
            "brca1_gallus_gallus",
            "brca1_homo_sapiens",
            "brca1_macaca_mulatta",
            "brca1_mus_musculus",
        ],
    ]

    # Expand the score matrix to include 'N' as a mismatch
    score_matrix["N"] = {base: 5 for base in "ACGTN"}

    with open("alignment.fasta", "w") as f:
        for i, seq_set in enumerate(seq_sets):
            seqs = {
                key: value
                for key, value in load_sequences(r"Projects\Project 3\data\brca1-testseqs.fasta").items()
                if key in seq_set
            }
            optimal, alignment = MSA(seqs, score_matrix, gap_penalty)
            print(f"Optimal score for {seq_set}: {optimal}")
            f.write(f"# Alignment {i+1}\n")
            for key, value in alignment.items():
                f.write(f">{key}\n{value}\n")

        # For the full length BRCA1 genes
        seqs = load_sequences(r"Projects\Project 3\data\brca1-full.fasta")
        optimal, alignment = MSA(seqs, score_matrix, gap_penalty)
        print(f"Optimal score for full length BRCA1 genes: {optimal}")
        f.write("# Full length BRCA1 genes alignment\n")
        for key, value in alignment.items():
            f.write(f">{key}\n{value}\n")


# presentation()
