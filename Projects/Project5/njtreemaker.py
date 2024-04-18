from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import numpy as np
from typing import Tuple


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


def select_neighbors(n_matrix, tree: Tree) -> tuple[Clade, Clade]:
	min_n_ij = 10000

	min_i = float("inf")
	min_j = float("inf")

	for i in range(1, len(n_matrix)):
		for j in range(1, len(n_matrix)):
			if i != j and n_matrix[i][j] < min_n_ij:
				min_n_ij = n_matrix[i][j]
				min_i = i
				min_j = j

	n_matrix[min_i][0], n_matrix[0][min_j]

	clade_i = next(tree.find_elements(name=n_matrix[min_i][0]))
	clade_j = next(tree.find_elements(name=n_matrix[0][min_j]))

	return clade_i, clade_j


def update_matrix(tree: Tree) -> list[list]:

	return matrix


def collapse_neighbors(tree: Tree, clades: tuple[Clade, Clade]) -> Tree:

	clade_i, clade_j = clades
	new_clade = Clade(name=f"k_{clade_i.name}{clade_j.name}")

	new_clade.clades.append(clade_i)
	new_clade.clades.append(clade_j)

	tree.root.clades.remove(clade_i)
	tree.root.clades.remove(clade_j)

	tree.root.clades.append(new_clade)

	return tree


def terminate(tree) -> Tree:

	return tree


def save_newick(tree: Tree, filename: str):
	Phylo.write(tree, filename, "newick")
	return


def main():
	matrix = load_distancematrix("Projects/Project5/example_slide4.phy")

	tree = build_root_tree(matrix)

	# print(n_matrix(matrix, tree))
	# print(select_neighbors(n_matrix(matrix, tree)))

	# print(tree)
	# Phylo.draw(tree)

	clades = select_neighbors(n_matrix(matrix, tree), tree)
	tree = collapse_neighbors(tree, clades)
	
	print(tree)


if __name__ == "__main__":
	main()
