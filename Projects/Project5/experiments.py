import timeit
import os
import subprocess

# Define the commands for njtreemaker, quicktree and rapidnj
# commands = [
#     "njtreemaker",
#     "quicktree",
#     "rapidnj",
# ]
commands = ["njtreemaker"]
# Define the path to the distance matrices
matrix_path = "Projects/Project5/unique_distance_matrices"

output_path = "Projects/Project5/unique_trees"

# Get a list of all distance matrices
matrices = os.listdir(matrix_path)

# For debugging, we'll just use the first matrix
matrix = matrices[0]

# Time each command
for command in commands:
    start_time = timeit.default_timer()
    match command:
        case "njtreemaker":
            subprocess.run(
                f" python3 Projects/Project5/njtreemaker.py {os.path.join(matrix_path, matrix)} {output_path}/njtreemaker_{os.path.basename(matrix)}.new",
                shell=True,
            )
        case "quicktree":
            subprocess.run(f"quicktree -in m", shell=True)
        case "rapidnj":
            subprocess.run(f"rapidNJ ", shell=True)

    subprocess.run(f"{command} ", shell=True)
    end_time = timeit.default_timer()
    print(f"{command} took {end_time - start_time} seconds to run on {matrix}")
