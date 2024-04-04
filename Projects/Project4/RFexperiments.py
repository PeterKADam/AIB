import RFDIST


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


def main(files: dict):
    # Initialize a dictionary to store the RF distances
    rf_matrix = {file: {file: 0 for file in files.values()} for file in files.values()}

    # Get the list of file names
    files = list(files.values())

    # Iterate over each pair of files
    for i in range(6):
        for j in range(i, 6):
            # Load the trees
            tree1 = RFDIST.loadtree(files[i])
            tree2 = RFDIST.loadtree(files[j])

            # Calculate the RF distance
            rf_distance = RFDIST.calcRFDist(tree1, tree2)

            # Store the RF distance in the matrix
            rf_matrix[files[i]][files[j]] = rf_distance
            rf_matrix[files[j]][files[i]] = rf_distance  # the matrix is symmetric

    # Print the RF distance matrix
    for file1, row in rf_matrix.items():
        print(f"{file1}:")
        for file2, distance in row.items():
            print(f"  {file2}: {distance}")


main(filesexp1)
