
import


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

def build_root_tree(matrix):
	return tree

def update_matrix(tree):
	return matrix

def add_node_to_tree(tree,something):
	return tree

def terminate(tree):
	return tree

def save_newick(tree):
	return
 