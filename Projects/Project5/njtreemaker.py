from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import numpy as np


def load_distancematrix(path):
    distance_matrix = []

    with open(path, "r") as file:
        for line_num, line in enumerate(file):
            line = line.strip()
            if line:
                elements = line.split()
                if line_num == 0:  # First line, indicating matrix size
                    size = int(elements[0])
                    labels = elements[1:]
                    distance_matrix.append([size] + labels)
                else:
                    label = elements[0]
                    distances = list(map(float, elements[1:]))
                    new_line = [label] + distances
                    distance_matrix[0].append(label)
                    distance_matrix.append(new_line)

    return distance_matrix


def build_root_tree(matrix: list[list]) -> Tree:

    root_clade = Clade()

    clades = [Clade(None, name) for name in matrix[0][1:]]

    for clade in clades:
        root_clade.clades.append(clade)

    tree = Tree(root=root_clade, rooted=False)

    return tree


def n_matrix(matrix, tree):

    n_matrix = [[0 for _ in range(len(matrix))] for _ in range(len(matrix))]

    n_matrix[0] = matrix[0]
    for i in range(1, len(matrix)):
        n_matrix[i][0] = matrix[i][0]

    for i in range(1, len(matrix)):
        for j in range(1, len(matrix)):

            d_ij = matrix[i][j]
            r_i = 1 / matrix[0][0] * sum(matrix[i][1:])
            r_j = 1 / matrix[0][0] * sum(matrix[j][1:])
            n_ij = d_ij - (r_i + r_j)
            n_matrix[i][j] = n_ij
    return n_matrix


def select_neighbors(n_matrix):
    min_ij = 10000

    min_i = float("inf")
    min_j = float("inf")

    for i in range(1, len(n_matrix)):
        for j in range(1, len(n_matrix)):
            if i != j and n_matrix[i][j] < min_ij:
                min_ij = n_matrix[i][j]
                min_i = i
                min_j = j

    return n_matrix[min_i][0], n_matrix[0][min_j]


def update_matrix(tree: Tree) -> list[list]:

    return matrix


def add_node_to_tree(tree: Tree, something) -> Tree:

    return tree


def terminate(tree) -> Tree:

    return tree


def save_newick(tree: Tree, filename: str):
    Phylo.write(tree, filename, "newick")
    return


matrix = load_distancematrix("Projects/Project5/example_slide4.phy")
tree = build_root_tree(matrix)

print(n_matrix(matrix, tree))
print(select_neighbors(n_matrix(matrix, tree)))

# print(tree)
# Phylo.draw(tree)
