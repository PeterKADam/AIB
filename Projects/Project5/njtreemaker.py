from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade


<<<<<<< HEAD
def load_distancematrix(path):  
	distance_matrix = []

    with open(path, "r") as file:
        for line_num, line in enumerate(file):
            line = line.strip()
            if line:
                elements = line.split()
                if line_num == 0:  # First line, indicating matrix size
                    size = int(elements[0])
                    labels = [""] + elements[1:]
                    distance_matrix.append([size] + labels)
                else:
                    label = elements[0]
                    distances = list(map(float, elements[1:]))
                    new_line = [label] + distances
                    distance_matrix.append(new_line)

    return distance_matrix

print(load_distancematrix("example_slide4.phy"))
=======
def load_distancematrix():
    return matrix
>>>>>>> 73e92150183b873873e431d6b8f3512fc5239026


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
