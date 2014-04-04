#!/usr/bin/env python2
#-*- coding: utf-8 -*-

"""
A custom made script to modify an OTU table for a specific project in a specific way.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./modify_otu_table ~/ILLUMITAG/views/clusters/danube/otus/uparse/taxonomy_fw/otu_table.csv

Or like this:
$ ./modify_otu_table ~/ILLUMITAG/views/clusters/danube/otus/uparse/taxonomy_silva/otu_table.csv
"""

###############################################################################
# The libraries we need #
import os, sys, pandas
# Check usage #
if len(sys.argv) != 2: sys.exit(sys.modules[__name__].__doc__)
# Get the shell arguments #
otu_path = os.path.expanduser(sys.argv[1])
# Check that the path is valid #
if not os.path.exists(otu_path): raise Exception("No file at %s." % otu_path)
# Read #
otus = pandas.io.parsers.read_csv(otu_path, sep='\t', index_col=0, encoding='utf-8')
# Define operations #
operations = {
    "II003F": ("II003Ftr2", "II003Ftr3"),
    "II009A": ("II009Abr2",),
    "II014A": ("II014Atr2", "II014Atr3"),
    "II015A": ("II015Abr2", "II015Atr2", "II015Atr3"),
    "II021A": ("II021Abr2",),
    "II023F": ("II023Ftr2", "II023Ftr3"),
    "II024A": ("II024Abr2", "II024Atr2", "II024Atr3"),
    "II032A": ("II032Atr2", "II032Atr3"),
    "II032F": ("II032Fbr2",),
    "II035F": ("II035Fbr2",),
    "II047A": ("II047Atr2", "II047Atr3"),
    "II052A": ("II052Abr2", "II052Atr3"),
    "II071A": ("II071Abr2",),
    "II071F": ("II071Fbr2",),
    "II078A": ("II078Atr2", "II078Atr3"),
    "II079A": ("II079Atr2", "II079Atr3"),
    "II084F": ("II084Ftr2", "II084Ftr3")
}
# Execute them #
for parent, children in operations.items():
    otus.loc[parent] = sum((otus.loc[child] for child in children), otus.loc[parent])
    for child in children: otus = otus.drop(child)
# Write #
otus.to_csv(otu_path, sep='\t')