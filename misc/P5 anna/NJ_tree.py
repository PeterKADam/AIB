from __future__ import division
import sys 
import decimal

def make_table(n, m):
    return [[float("+inf") for i in range(m)] for j in range(n)]

def matrix_N(distances, seqs_left, tree, seen, seqs):

	N = make_table(len(distances),len(distances[0]))
	minimum = {}

	for i in range(len(distances)):

		ri = (float(1/(seqs_left-2)) * sum( [ float(distances[i][x]) if distances[i][x] != float('-inf') else 0 for x in range(len(distances[i])) ] ))

		for j in range(len(distances)):

			if i != j and N[j][i] == float('+inf') and distances[i][j] != float('-inf'):

				rj = (float(1/(seqs_left-2)) * sum( [ float(distances[j][x]) if distances[j][x] != float('-inf') else 0 for x in range(len(distances[j])) ]))

				N[i][j] = float(distances[i][j]) - (ri + rj)

				if i in list(minimum.keys()):

					if N[i][j] < minimum[i][0]:

						minimum[N[i][j]] = [i,j]
						weight_i = 0.5 * (distances[i][j] + ri - rj)
						weight_j = 0.5 * (distances[i][j] + rj - ri)

				else:

					weight_i = 0.5 * (float(distances[i][j]) + ri - rj)
					weight_j = 0.5 * (float(distances[i][j]) + rj - ri)
					minimum[N[i][j]] = [[i,j],[weight_i, weight_j]]

	minimum_all = [v for k,v in minimum.items() if k == min(list(minimum.keys()))][0]
	get_tree = append_tree(tree, minimum_all, seqs_left, distances, seen, seqs)

	tree = get_tree[0]
	seqs_left = get_tree[1]

	i,j = [int(x) for x in minimum_all[0]]

	for m in range(len(distances)):

		if m != i and m != j:

			if distances[m][i] != float('-inf'):

				distances[m].append( 0.5 * (float(distances[m][i])+float(distances[m][j]) - float(distances[i][j])) )

			else:

				distances[m].append(float('-inf'))
		else:

			distances[m].append(float('-inf'))

	distances.append([ 0.5  * (float(distances[i][m])+float(distances[j][m]) - float(distances[i][j])) if distances[i][m] != float('-inf') else float('-inf') for m in range(len(distances)+1) ])
	distances[i] = distances[j] = [float('-inf') for x in range(len(distances))]


	return N,minimum_all,tree,distances,seqs_left,seen

def append_tree(tree, minimum_all, seqs_left, distances, seen, seqs):

	group = ()

	leaves = minimum_all[0]
	print(leaves)
	weights = minimum_all[1]
	print(weights)
	to_remove = []

	for i in range(len(leaves)):

		index = leaves[i]
		
		if leaves[i] in tree.keys():

			leaf = (":"+str(weights[i]),)
			group = group + (tree[leaves[i]],)+leaf
			to_remove.append(leaves[i])

		else:

			group = group + (str(seqs[leaves[i]])+":"+str(weights[i]),)
			seqs_left-=1
			seen.append(leaves[i])

	if len(list(tree.keys())) > 0:
	
		new_node = max(list(tree.keys())) + 1

	else:

		new_node = len(distances)

	tree[new_node] = group

	for i in to_remove:
		tree.pop(i)

	return tree,seqs_left,seen

in_file = sys.argv[1]
out_file_2 = sys.argv[2]

distances = []
seqs = []

with open(in_file, 'r') as f:
	lines = f.read().rstrip('\n').split('\n')
	for i in range(len(lines)):
		
		if i != 0:
			elems = lines[i].split(" ")
			distances.append(elems[1:])
			seqs.append(elems[0])			

f.close()

all_seqs = len(distances)
tree = {}
nodes = 0
new_node = len(distances)
node_index_dict	= {}

seen = []

seqs_left = all_seqs

while seqs_left > 3:

	test = matrix_N(distances, seqs_left, tree, seen, seqs)
	N = test[0]
	minimum = test[1]
	tree = test[2]
	distances = test[3]
	seqs_left = test[4]
	seen = test[5]
	print(seqs_left)

new_node = max(list(tree.keys())) + 1
to_finish = []

for i in range(all_seqs):

	if not i in seen:

		to_finish.append(i)
print(to_finish)

i = to_finish[0]
j = to_finish[1]
m = to_finish[2]

wi = (float(distances[i][j]) + float(distances[i][m]) - float(distances[j][m])) / 2
wj = (float(distances[i][j]) + float(distances[j][m]) - float(distances[i][m])) / 2
wm = (float(distances[i][m]) + float(distances[j][m]) - float(distances[i][j])) / 2

tree[new_node]=(seqs[i]+":"+str(wi), seqs[j]+":"+str(wj), seqs[m]+":"+str(wm))

tree_newick = ()

for i in list(tree.keys()):

	tree_newick = tree_newick + (tree[i],)

print(tree_newick)

out_file = out_file_2+".tmp"

with open(out_file, "w") as o:
	o.write(str(tree_newick))
	o.write(';')

import subprocess

args1 = 'tr -d [" \'"] < ' + out_file + '|awk -F ",:" ' + "'{OFS=" + ' ":" ' + ";$1=$1; print $0}' > " + out_file_2
subprocess.call(args1, shell=True)

import os
os.remove(out_file)

		


