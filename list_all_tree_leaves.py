#!/usr/bin/env python3

import sys

import dendropy

tree = dendropy.Tree.get(path=sys.argv[1], schema='newick')

for leaf in tree.taxon_namespace:
    print(leaf)
