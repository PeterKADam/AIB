import csv
import timeit
import subprocess
import os

# Define the commands for njtreemakernumpy, quicktree and rapidnj
commands = ["njtreemakernumpy", "quicktree", "rapidnj", "njtreemakernumpy"]
# commands = ["njtreemakernumpy"]
# Define the path to the distance matrices
matrix_path = "Projects/Project5/unique_distance_matrices"

output_path = "Projects/Project5/unique_trees"

# Get a list of all distance matrices
matrices = os.listdir(matrix_path)

# matrix = [matrices[-1]]
# Open the TSV file in write mode
with open("runtimes.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    # Write the header
    writer.writerow(
        [
            "Matrix",
            "QuickTree Time",
            "RapidNJ Time",
            "njtreemakernumpy",
            "Speed-up QT",
            "Speed-up RNJ",
            "RF QT-NJTM",
            "RF RNJ-NJTM",
            "RF RNJ-QT",
        ]
    )

    # Time each command
    for matrix in matrices[:-1]:
        times = {}
        trees = {}
        for command in ["quicktree", "rapidnj", "njtreemakernumpy"]:
            start_time = timeit.default_timer()
            if command == "njtreemakernumpy":
                tree_file = (
                    f"{output_path}/njtreemakernumpy_{os.path.basename(matrix)}.new"
                )
                subprocess.run(
                    f"python3 Projects/Project5/njtreemakernumpy.py {os.path.join(matrix_path, matrix)} {tree_file}",
                    shell=True,
                )
            elif command == "quicktree":
                tree_file = f"{output_path}/quicktree_{os.path.basename(matrix)}.new"
                subprocess.run(
                    f"quicktree -in m -out t {os.path.join(matrix_path, matrix)} > {tree_file}",
                    shell=True,
                )
            elif command == "rapidnj":
                tree_file = f"{output_path}/rapidnj_{os.path.basename(matrix)}.new"
                subprocess.run(
                    f"rapidnj {os.path.join(matrix_path, matrix)} -i pd -o t -x {tree_file}",
                    shell=True,
                )
            end_time = timeit.default_timer()
            times[command] = end_time - start_time
            trees[command] = tree_file

        # Compute RF distances
        rf_distances = {}
        for pair in [
            ("quicktree", "njtreemakernumpy"),
            ("rapidnj", "njtreemakernumpy"),
            ("rapidnj", "quicktree"),
        ]:
            output = subprocess.check_output(
                f"python3 Projects/Project4/RFDIST.py {trees[pair[0]]} {trees[pair[1]]}",
                shell=True,
            )
            rf_distances[pair] = float(output)

        # Write the command, matrix, and runtime to the TSV file
        writer.writerow(
            [
                matrix,
                times["quicktree"],
                times["rapidnj"],
                times["njtreemakernumpy"],
                times["quicktree"] / times["njtreemakernumpy"],
                times["rapidnj"] / times["njtreemakernumpy"],
                rf_distances[("quicktree", "njtreemakernumpy")],
                rf_distances[("rapidnj", "njtreemakernumpy")],
                rf_distances[("rapidnj", "quicktree")],
            ]
        )
