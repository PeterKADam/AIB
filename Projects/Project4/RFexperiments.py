import RFDIST

prepath = "Projects/Project4/"

filesexp1 = {
    "qt_muscle": "QT_Exp1_Muscle.stockholm.new",
    "qt_kalign": "QT_Exp1_Kalign.stockholm.new",
    "qt_omega": "QT_Exp1_Omega.stockholm.new",
    "RNJ_muscle": "RNJ_Exp1_Muscle.stockholm.new",
    "RNJ_kalign": "RNJ_Exp1_Kalign.stockholm.new",
    "RNJ_omega": "RNJ_Exp1_Omega.stockholm.new",
}

filesexp2 = {
    "qt_muscle": "QT_Exp2_Muscle.stockholm.new",
    "qt_kalign": "QT_Exp2_Kalign.stockholm.new",
    "qt_omega": "QT_Exp2_Omega.stockholm.new",
    "RNJ_muscle": "RNJ_Exp2_Muscle.stockholm.new",
    "RNJ_kalign": "RNJ_Exp2_Kalign.stockholm.new",
    "RNJ_omega": "RNJ_Exp2_Omega.stockholm.new",
}


def main(files: dict, prepath: str):
    # Initialize a dictionary to store the RF distances
    rf_matrix = {file: {file: 0 for file in files.values()} for file in files.values()}

    # Get the list of file names
    files = list(files.values())

    # Iterate over each pair of files
    for i in range(6):
        for j in range(i, 6):
            # Load the trees, prepending the prepath to the file name
            tree1 = RFDIST.loadtree(prepath + files[i])
            tree2 = RFDIST.loadtree(prepath + files[j])

            # Find the splits in each tree
            splits1 = RFDIST.findsplits(tree1)
            splits2 = RFDIST.findsplits(tree2)

            # Calculate the RF distance
            rf_distance = RFDIST.calcRFDist(splits1, splits2)

            # Store the RF distance in the matrix
            rf_matrix[files[i]][files[j]] = rf_distance
            rf_matrix[files[j]][files[i]] = rf_distance

    # Print the RF distance matrix
    print(" " * 10, end=" ")  # padding for alignment
    for file in files:
        print(
            f"{file[:10]:>10}", end=" "
        )  # print the first 10 characters of each file name
    print()

    for file1, row in rf_matrix.items():
        print(
            f"{file1[:10]:<10}", end=""
        )  # print the first 10 characters of each file name
        for file2, distance in row.items():
            print(f"{distance:>10}", end="")
        print()


# main(filesexp1, prepath)

# main(filesexp2, prepath)


def main2(files1: dict, files2: dict, prepath: str):
    # Initialize a dictionary to store the RF distances
    rf_distances = {}

    # Iterate over each method
    for method in files1.keys():
        # Load the trees, prepending the prepath to the file name
        tree1 = RFDIST.loadtree(prepath + files1[method])
        tree2 = RFDIST.loadtree(prepath + files2[method])

        # Find the splits in each tree
        splits1 = RFDIST.findsplits(tree1)
        splits2 = RFDIST.findsplits(tree2)

        # Calculate the RF distance
        rf_distance = RFDIST.calcRFDist(splits1, splits2)

        # Store the RF distance
        rf_distances[method] = rf_distance

    # Print the RF distances
    for method, distance in rf_distances.items():
        print(f"{method[:10]:<10}: {distance}")


main2(filesexp1, filesexp2, prepath)
