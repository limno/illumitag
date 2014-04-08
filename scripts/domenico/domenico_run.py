#!/usr/bin/env python2

"""
A script to run the clustering analyses.
"""

# Future #
from __future__ import division

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag

###############################################################################
# Create the cluster #
# cluster = illumitag.clustering.favorites.domenico
cluster = illumitag.clustering.favorites.danube

# Reporting #
cluster.reporter.total_raw_count
cluster.reporter.total_seq_count
cluster.reporter.otu_count_raw
cluster.reporter.otu_count

# New graph #
cluster.otu_uparse.taxonomy_silva.graphs[-3].plot()

# Special values #
targets = ['LD12', 'acI-B1', 'acI-A7', 'acI-C2']
table = cluster.otu_uparse.taxonomy_fw.comp_tips.taxa_table
cols = table.columns.isin(targets)
ratio = table[table.columns[cols]].sum(axis=1) / table[table.columns[~cols]].sum(axis=1)
ratio.sort()