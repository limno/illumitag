#!/usr/bin/env python2

"""
A script to run the clustering analyses.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag
from tqdm import tqdm

###############################################################################
# Create the cluster #
cluster = illumitag.clustering.Cluster(illumitag.runs[0][0].samples.children, 'test')

# Run it #
illumitag.runs[0][0].create_samples()
for s in tqdm(cluster): s.process()
cluster.combine_reads()
cluster.set_size(400)
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()
cluster.otu_uparse.taxonomy_silva.make_otu_table()
cluster.otu_uparse.taxonomy_silva.make_plots()
cluster.otu_uparse.taxonomy_silva.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_uparse.taxonomy_silva.comp_phyla.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_plots()