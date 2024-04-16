class alignment_matrix:
    def __init__(self, sequence1, sequence2, score_matrix, gap_cost):
        self.sequence1 = sequence1.upper()
        self.sequence2 = sequence2.upper()
        self.score_matrix = score_matrix
        self.gap_cost = gap_cost
        self.matrix = self.empty_matrix()
        self.alignments = []
        self.completealignment()

    def get_max(self):
        return self.matrix[-1][-1]

    def completealignment(self):
        self.prepare_matrixes()
        self.fill_matrixes()
        # self.get_alignments()

    def empty_matrix(self):
        return [[None] * (len(self.sequence1) + 1) for _ in range(len(self.sequence2) + 1)]

    def print_dp_matrix(self, matrix_object):
        if len(matrix_object) > 20:
            return print("Matrix too large to print")

        max_len = max(len(str(cell)) for row in matrix_object for cell in row)
        fmt = "{{:>{}}}".format(max_len + 1)
        row_fmt = fmt * (len(matrix_object[0]) + 1) + "\n"
        mat_fmt = row_fmt * (len(matrix_object) + 1)
        seq1, seq2 = self.sequence1, self.sequence2
        seq1 = " " + seq1
        seq2 = " " + seq2
        lst = [" "] + list(seq1)
        for i in range(len(seq2)):
            lst.extend([seq2[i]] + list(map(repr, matrix_object[i])))
        print(mat_fmt.format(*lst))

    def niceprint_alignments(self):
        if len(self.alignments) == 0:
            print("No alignments found")
        elif len(self.alignments) <= 5:
            for i, a in enumerate(self.alignments, start=1):
                print(f"Pair {i}:\n\t{a[0]}\n\t{a[1]}\n")
        else:
            print(
                f"Too many alignments to print, there are: {len(self.alignments)} alignments\n"
                f"the start of one of the alignments is:\n\t{self.alignments[0][0][:100]}\n\t{self.alignments[0][1][:100]}\n"
            )

    def get_alignments(self):
        alignments = self.traceback_recursive(len(self.sequence2), len(self.sequence1))
        # Reverse the completed strings
        self.alignments = [(a1[::-1], a2[::-1]) for a1, a2 in alignments]
        return self.alignments

    def get_alignment(self):
        alignments = self.traceback_nonrecursive(len(self.sequence2), len(self.sequence1))
        # Reverse the completed strings
        self.alignments = [(a1[::-1], a2[::-1]) for a1, a2 in alignments]
        return self.alignments

    def print_finish_matrix(self):

        self.print_dp_matrix(self.matrix)
        print(f"maximum is : {self.matrix[-1][-1]}")
        self.get_alignments()
        self.niceprint_alignments()

    def traceback_recursive(self, row, col):
        if row == 0 and col == 0:
            return [("", "")]
        elif row == 0:
            return [(col * "-", self.sequence1[:col])]
        elif col == 0:
            return [(self.sequence2[:row], row * "-")]

        arrows = self.get_traceback_arrows(
            row,
            col,
        )

        alignments = []

        for arrow in arrows:
            n, direction = arrow.split("_")
            n = int(n)
            if direction == "diagonal":
                sub_alignments = self.traceback_recursive(row - n, col - n)
                alignments.extend(
                    [(self.sequence2[row - n] + a1, self.sequence1[col - n] + a2) for a1, a2 in sub_alignments]
                )
            elif direction == "up":
                sub_alignments = self.traceback_recursive(row - n, col)
                alignments.extend([(self.sequence2[row - n] + a1, (n * "-") + a2) for a1, a2 in sub_alignments])
            elif direction == "left":
                sub_alignments = self.traceback_recursive(row, col - n)
                alignments.extend([((n * "-") + a1, self.sequence1[col - n] + a2) for a1, a2 in sub_alignments])

        return alignments

    def traceback_nonrecursive(self, row, col):
        if row == 0 and col == 0:
            return [("", "")]
        elif row == 0:
            return [(col * "-", self.sequence1[:col])]
        elif col == 0:
            return [(self.sequence2[:row], row * "-")]

        arrows = [
            self.get_traceback_arrows(
                row,
                col,
            )[0]
        ]  # now this is software engineering

        alignments = []

        for arrow in arrows:
            n, direction = arrow.split("_")
            n = int(n)
            if direction == "diagonal":
                sub_alignments = self.traceback_nonrecursive(row - n, col - n)
                alignments.extend(
                    [(self.sequence2[row - n] + a1, self.sequence1[col - n] + a2) for a1, a2 in sub_alignments]
                )
            elif direction == "up":
                sub_alignments = self.traceback_nonrecursive(row - n, col)
                alignments.extend([(self.sequence2[row - n] + a1, (n * "-") + a2) for a1, a2 in sub_alignments])
            elif direction == "left":
                sub_alignments = self.traceback_nonrecursive(row, col - n)
                alignments.extend([((n * "-") + a1, self.sequence1[col - n] + a2) for a1, a2 in sub_alignments])

        return alignments

    def get_traceback_arrows(self, row, col):
        score_diagonal = self.matrix[row - 1][col - 1]

        score_current = self.matrix[row][col]

        arrows = []

        match_score = self.score_matrix[self.sequence2[row - 1]][self.sequence1[col - 1]]

        if score_current == score_diagonal + match_score:
            arrows.append("1_diagonal")

        for n in range(1, row, 1):
            if score_current == self.matrix[row - n][col] + self.gap_open + (self.gap_cost * n):
                arrows.append(f"{n}_up")
                break

        for n in range(1, col, 1):
            if score_current == self.matrix[row][col - n] + self.gap_open + (self.gap_cost * n):
                arrows.append(f"{n}_left")
                break

        return arrows


class LinearGlobalAlignment(alignment_matrix):
    def __init__(self, sequence1, sequence2, score_matrix, gap_cost):
        self.gap_open = 0
        super().__init__(sequence1, sequence2, score_matrix, gap_cost)

    def prepare_matrixes(self):
        for col in range(len(self.sequence1) + 1):
            self.matrix[0][col] = col * self.gap_cost

        for row in range(len(self.sequence2) + 1):
            self.matrix[row][0] = row * self.gap_cost

    def fill_matrixes(self):

        for col in range(1, len(self.sequence1) + 1):  # j
            for row in range(1, len(self.sequence2) + 1):  # i

                self.matrix[row][col] = min(
                    # C(i-1, j-1) + match
                    self.matrix[row - 1][col - 1] + self.score_matrix[self.sequence1[col - 1]][self.sequence2[row - 1]],
                    # c(i, j-1) + gapcost
                    self.matrix[row][col - 1] + self.gap_cost,
                    # c(i-1, j) + gapcost
                    self.matrix[row - 1][col] + self.gap_cost,
                )
