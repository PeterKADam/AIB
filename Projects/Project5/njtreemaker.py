from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import sys
import time
from functools import cache
from functools import lru_cache


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


def save_newick(tree: Tree, filename: str):
	Phylo.write(tree, filename, "newick")


def build_root_tree(matrix: list[list]) -> Tree:

	root_clade = Clade()

	for clade in [Clade(None, name) for name in matrix[0][1:]]:
		root_clade.clades.append(clade)

	return Tree(root=root_clade, rooted=False)



def r_x(matrix, x, cache={}):
	if x in cache:
		return cache[x]
	result = 1 / (matrix[0][0] - 2) * sum(matrix[x][1:])
	cache[x] = result
	return result

def n_matrix(matrix):
	
	def calc_nij(matrix, i, j, r_cache):
		d_ij = matrix[i][j]
		r_i = r_x(matrix, i, r_cache)
		r_j = r_x(matrix, j, r_cache)
		n_ij = d_ij - (r_i + r_j)
		return n_ij
	
	r_cache = {}
	n_matrix = [[0 for _ in range(len(matrix))] for _ in range(len(matrix))]

	n_matrix[0] = matrix[0]

	for i in range(1, len(matrix)):
		n_matrix[i][0] = matrix[i][0]

	for i in range(1, len(matrix)):
		for j in range(1, len(matrix)):
			if i != j:
				n_matrix[i][j] = calc_nij(matrix, i, j, r_cache)
			else:
				n_matrix[i][j] = float("inf")

	return n_matrix


def select_neighbors(n_matrix) -> tuple[Clade, Clade]:
	min_n_ij = float("inf")

	min_i = float("inf")
	min_j = float("inf")

	for i in range(1, len(n_matrix)):
		for j in range(1, len(n_matrix)):
			if i != j and n_matrix[i][j] < min_n_ij:
				min_n_ij = n_matrix[i][j]
				min_i = i
				min_j = j

	n_matrix[min_i][0], n_matrix[0][min_j]

	return min_i, min_j


def get_neighbor_clades(
	tree: Tree, n_matrix: list[list], min_i: int, min_j: int
) -> tuple[Clade, Clade]:

	clade_i = next(tree.find_elements(name=n_matrix[min_i][0]))
	clade_j = next(tree.find_elements(name=n_matrix[0][min_j]))

	return clade_i, clade_j


def update_d_matrix(
	matrix: list[list], clades: tuple[Clade, Clade], i: int, j: int
) -> list[list]:
	# Calculate k_ij for all m in S except {i, j}
	k_ij = [
		0.5 * (matrix[i][m] + matrix[j][m] - matrix[i][j])
		for m in range(1, len(matrix))  # Start from 1 to skip the header
		if m not in {i, j}
	]

	k_ij.append(0.0)

	# Remove the rows and columns corresponding to i and j
	matrix = [row for k, row in enumerate(matrix) if k not in {i, j}]
	matrix = [[col for l, col in enumerate(row) if l not in {i, j}] for row in matrix]

	clade_i, clade_j = clades
	matrix.append([f"{clade_i.name}{clade_j.name}"] + k_ij)

	for idx, row in enumerate(matrix[1:-1], start=1):
		row.append(k_ij[idx - 1])

	matrix[0].append(f"{clade_i.name}{clade_j.name}")

	matrix[0][0] -= 1  # |S|

	return matrix


def collapse_neighbors(
	tree: Tree, clades: tuple[Clade, Clade], n_matrix, i, j, matrix
) -> Tree:

	clade_i, clade_j = clades
	new_clade = Clade(name=f"{clade_i.name}{clade_j.name}")

	for clade in clades:
		new_clade.clades.append(clade)
		tree.root.clades.remove(clade)

	tree.root.clades.append(new_clade)

	gamma_ki = 0.5 * (matrix[i][j] + r_x(matrix, i) - r_x(matrix, j))
	gamma_kj = matrix[i][j] - gamma_ki

	clade_i.branch_length = gamma_ki
	clade_j.branch_length = gamma_kj

	return tree


def terminate(tree: Tree, matrix) -> Tree:

	v = tree.root
	v.name = "v"

	gamma_vi = (matrix[1][2] + matrix[1][3] - matrix[2][3]) / 2
	gamma_vj = (matrix[1][2] + matrix[2][3] - matrix[1][3]) / 2
	gamma_vm = (matrix[1][3] + matrix[2][3] - matrix[1][2]) / 2

	for clade in v.clades:
		if clade.name == matrix[1][0]:
			clade.branch_length = gamma_vi
		if clade.name == matrix[2][0]:
			clade.branch_length = gamma_vj
		if clade.name == matrix[3][0]:
			clade.branch_length = gamma_vm

	return tree


def main(input, output):
	matrix = load_distancematrix(input)

	tree = build_root_tree(matrix)

	while matrix[0][0] > 3:
		S = matrix[0][0]
		t1 = time.perf_counter()
		# print(f"S: {matrix[0][0]}")
		n = n_matrix(matrix)
		neighbors_index = select_neighbors(n)
		neighbors = get_neighbor_clades(tree, n, *neighbors_index)
		tree = collapse_neighbors(tree, neighbors, n, *neighbors_index, matrix)
		matrix = update_d_matrix(matrix, neighbors, *neighbors_index)
		t2 = time.perf_counter()
		print(f"S:{S}\tTime: {t2-t1:.2f} sec")

	tree = terminate(tree, matrix)
	# Phylo.draw(tree)

	Phylo.write(tree, output, "newick")


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
