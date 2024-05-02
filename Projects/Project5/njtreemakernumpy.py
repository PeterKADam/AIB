import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import sys


def load_distancematrix(path):
    with open(path, "r") as file:
        lines = file.readlines()
        headers = [line.strip().split()[0] for line in lines[1:]]  # Skip the first line
        matrix = np.array([list(map(float, line.strip().split()[1:])) for line in lines[1:]])  # Skip the first line
    return headers, matrix


def build_root_tree(headers: list) -> Tree:

    root_clade = Clade()

    for clade in [Clade(None, name) for name in headers]:
        root_clade.clades.append(clade)

    return Tree(root=root_clade, rooted=False)


def r_x(matrix, x, cache={}):
    if x in cache:
        return cache[x]
    result = 1 / (matrix.shape[0] - 2) * np.sum(matrix[x, :])
    cache[x] = result
    return result


def n_matrix(matrix):
    def calc_nij(matrix, i, j, r_cache):
        d_ij = matrix[i, j]
        r_i = r_x(matrix, i, r_cache)
        r_j = r_x(matrix, j, r_cache)
        n_ij = d_ij - (r_i + r_j)
        return n_ij

    r_cache = {}
    n_matrix = np.zeros(matrix.shape)

    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            n_matrix[i, j] = calc_nij(matrix, i, j, r_cache)
            n_matrix[j, i] = n_matrix[i, j]  # Use symmetry to avoid redundant calculations

    return n_matrix


def n_matrix_np(matrix):
    def r_x(matrix, r_cache):
        if r_cache.get(matrix.shape[0]) is None:
            r_cache[matrix.shape[0]] = matrix.sum(axis=1) / (2 * (matrix.shape[0] - 2))
        return r_cache[matrix.shape[0]]

    r_cache = {}
    r_values = r_x(matrix, r_cache)
    n_matrix = matrix - r_values[:, None] - r_values[None, :]
    np.fill_diagonal(n_matrix, 0)
    return n_matrix


def select_neighbors(n_matrix):
    min_n_ij = np.inf
    min_i = min_j = np.inf

    for i in range(n_matrix.shape[0]):
        for j in range(i + 1, n_matrix.shape[1]):
            if n_matrix[i, j] < min_n_ij:
                min_n_ij = n_matrix[i, j]
                min_i = i
                min_j = j

    return min_i, min_j


def get_neighbor_clades(tree: Tree, headers: list, min_i: int, min_j: int) -> tuple[Clade, Clade]:
    clade_i = next(tree.find_elements(name=headers[min_i]))
    clade_j = next(tree.find_elements(name=headers[min_j]))

    return clade_i, clade_j


def update_d_matrix(matrix, headers, i, j):
    # Calculate k_ij for all m in S except {i, j}
    k_ij = 0.5 * (matrix[i, :] + matrix[j, :] - matrix[i, j])
    k_ij = np.delete(k_ij, [i, j])

    # Remove the rows and columns corresponding to i and j
    matrix = np.delete(matrix, [i, j], axis=0)
    matrix = np.delete(matrix, [i, j], axis=1)

    # Append k_ij to the matrix
    matrix = np.append(matrix, k_ij.reshape(-1, 1), axis=1)
    matrix = np.append(matrix, np.append(k_ij, 0).reshape(1, -1), axis=0)

    # Update headers
    merged_clade_name = f"{headers[i]}{headers[j]}"
    headers = [headers[k] for k in range(len(headers)) if k != i and k != j]
    headers.append(merged_clade_name)

    return matrix, headers


def update_d_matrix_np(matrix, headers, i, j):
    # Create a mask for all indices except i and j
    mask = np.ones(matrix.shape[0], dtype=bool)
    mask[[i, j]] = False

    # Calculate k_ij for the necessary indices
    k_ij = 0.5 * (matrix[i, mask] + matrix[j, mask] - matrix[i, j])

    # Remove the rows and columns corresponding to i and j
    matrix = matrix[mask][:, mask]

    # Create a new row and column for k_ij
    new_row = np.append(k_ij, 0).reshape(1, -1)
    new_col = np.append(k_ij, [0]).reshape(-1, 1)

    # Append the new row and column to the matrix
    matrix = np.hstack((matrix, new_col))
    matrix = np.vstack((matrix, new_row))

    # Update headers
    merged_clade_name = f"{headers[i]}{headers[j]}"
    headers = [headers[k] for k in range(len(headers)) if k != i and k != j]
    headers.append(merged_clade_name)

    return matrix, headers


def collapse_neighbors(tree: Tree, headers: list, min_i: int, min_j: int, n_matrix, matrix) -> Tree:
    clade_i, clade_j = get_neighbor_clades(tree, headers, min_i, min_j)
    new_clade = Clade(name=f"{clade_i.name}{clade_j.name}")

    for clade in [clade_i, clade_j]:
        new_clade.clades.append(clade)
        tree.root.clades.remove(clade)

    tree.root.clades.append(new_clade)

    gamma_ki = 0.5 * (matrix[min_i][min_j] + r_x(matrix, min_i, {}) - r_x(matrix, min_j, {}))
    gamma_kj = matrix[min_i][min_j] - gamma_ki

    clade_i.branch_length = gamma_ki
    clade_j.branch_length = gamma_kj

    return tree


def terminate(tree: Tree, headers: list, matrix) -> Tree:
    v = tree.root
    v.name = "v"

    gamma_vi = (matrix[0][1] + matrix[0][2] - matrix[1][2]) / 2
    gamma_vj = (matrix[0][1] + matrix[1][2] - matrix[0][2]) / 2
    gamma_vm = (matrix[0][2] + matrix[1][2] - matrix[0][1]) / 2

    for clade in v.clades:
        if clade.name == headers[0]:
            clade.branch_length = gamma_vi
        if clade.name == headers[1]:
            clade.branch_length = gamma_vj
        if clade.name == headers[2]:
            clade.branch_length = gamma_vm

    return tree


import time


def main(input, output):
    headers, matrix = load_distancematrix(input)

    tree = build_root_tree(headers)

    while matrix.shape[0] > 3:

        S = matrix.shape[0]
        t1 = time.perf_counter()
        n = n_matrix_np(matrix)
        t2 = time.perf_counter()
        t3 = time.perf_counter()
        n2 = n_matrix(matrix)
        t4 = time.perf_counter()
        min_i, min_j = select_neighbors(n)
        tree = collapse_neighbors(tree, headers, min_i, min_j, n, matrix)
        matrix, headers = update_d_matrix(matrix, headers, min_i, min_j)

        print(f"S:{S}\tTime: {t2-t1:.3f} | {t4-t3:.3f} sec")

    tree = terminate(tree, headers, matrix)

    Phylo.write(tree, output, "newick")


# main("Projects/Project5/example_slide4.phy", "output.newick")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

# The rest of the code remains the same, but you'll need to adjust the functions that use the headers to use the indices instead.
