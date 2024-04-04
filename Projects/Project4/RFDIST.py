from Bio import Phylo
from Bio.Phylo import BaseTree

from typing import List, Tuple, Set
import sys


def loadtree(path: str) -> BaseTree.Tree:
    return Phylo.read(path, "newick")


def findsplits(tree: BaseTree.Tree) -> List[Tuple[Set[str], Set[str]]]:
    def get_leaf_names(clade: BaseTree.Clade) -> Set[str]:
        if clade.is_terminal():  # base case: if the clade is a leaf
            return {clade.name}
        else:  # recursive case: if the clade is an internal node
            leaf_names = set()
            for child in clade.clades:
                leaf_names.update(get_leaf_names(child))
        return leaf_names

    all_leaves = get_leaf_names(tree.root)  # get all leaf names in the tree
    splits = []
    for clade in tree.find_clades(
        terminal=False
    ):  # iterate over all clades, not just direct children of root
        if not clade.is_terminal():
            leaf_names = get_leaf_names(clade)  # get leaf names in this clade
            other_leaves = all_leaves - leaf_names  # get leaf names not in this clade
            if (
                leaf_names and other_leaves
            ):  # only add the split if neither side is empty
                splits.append(
                    (leaf_names, other_leaves)
                )  # add pair of splits to the list
    return splits


def calcRFDist(
    splits1: List[Tuple[Set[str], Set[str]]], splits2: List[Tuple[Set[str], Set[str]]]
) -> int:
    shared = 0
    for split1 in splits1:
        for split2 in splits2:
            if split1 == split2:
                shared += 1
    return len(splits1) + len(splits2) - 2 * shared


# t1 = loadtree(sys.argv[1])
# t2 = loadtree(sys.argv[2])

# print(calcRFDist(findsplits(t1), findsplits(t2)))
