from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade


def load_distancematrix():
    return matrix


def build_root_tree(matrix: list[list]) -> Tree:

    clades = [Clade(None, name) for name in matrix[0][1:]]

    tree = Tree(clades=clades, rooted=False)

    return tree


def update_matrix(tree: Tree) -> list[list]:
    return matrix


def add_node_to_tree(tree: Tree, something) -> Tree:
    return tree


def terminate(tree) -> Tree:
    return tree


def save_newick(tree: BaseTree.Tree, filename: str):
    Phylo.write(tree, filename, "newick")
    return
