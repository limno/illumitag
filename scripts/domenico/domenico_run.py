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
import illumitag, pandas
from collections import defaultdict

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
ratio = table[table.columns[cols]].sum(axis=1) / table.sum(axis=1)
ratio.sort()

### Algae analysis ###
# Make new unfiltered taxa table #
otus = cluster.otu_uparse.taxonomy_silva.cluster_counts_table
taxa = defaultdict(lambda: defaultdict(int))
for sample_name, column in otus.iterrows():
    for otu_name, count in column.iteritems():
        species = cluster.otu_uparse.taxonomy_silva.assignments[otu_name]
        taxa[species][sample_name] += count
taxa = pandas.DataFrame(taxa).fillna(0).astype(int)
# Reads fraction per sample that are plastids #
cols = pandas.Series([len(col) > 2 and col[2] == 'Plastid' for col in taxa.columns])
ratio = taxa[taxa.columns[cols]].sum(axis=1) / taxa.sum(axis=1)
ratio.to_csv(sys.stdout, sep='\t', float_format='%.5g')
# Breakdown #
samples = [s.short_name for s in cluster.samples if s.info['Filter_fraction'] == '3.0']
plastids = taxa[taxa.columns[cols]].loc[samples]
for i,row in plastids.iterrows(): print '%s\t%s\t%s' % (row.name, row.idxmax(), row.max())