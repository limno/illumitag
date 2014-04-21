#!/usr/bin/env python2

"""
A script to run the clustering analyses.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag

###############################################################################
# Create the cluster #
cluster = illumitag.clustering.favorites.test

# Run it with UPARSE #
illumitag.runs[0][0].create_samples()
cluster.process_samples()
cluster.combine_reads()
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()
cluster.otu_uparse.taxonomy_silva.make_otu_table()
cluster.otu_uparse.taxonomy_silva.make_otu_table_norm()
cluster.otu_uparse.taxonomy_silva.make_plots()
cluster.otu_uparse.taxonomy_silva.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_uparse.taxonomy_silva.comp_phyla.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_plots()

# Run it with UCLUST #
cluster.otu_uclust.run()
cluster.otu_uclust.taxonomy_silva.assign()
cluster.otu_uclust.taxonomy_silva.make_otu_table()
cluster.otu_uclust.taxonomy_silva.make_otu_table_norm()
cluster.otu_uclust.taxonomy_silva.make_plots()
cluster.otu_uclust.taxonomy_silva.stats.nmds.run()
cluster.otu_uclust.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_uclust.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_uclust.taxonomy_silva.comp_phyla.stats.nmds.run()
cluster.otu_uclust.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_uclust.taxonomy_silva.comp_tips.make_plots()

# Run it with CDHIT #
cluster.otu_cdhit.run()
cluster.otu_cdhit.taxonomy_silva.assign()
cluster.otu_cdhit.taxonomy_silva.make_otu_table()
cluster.otu_cdhit.taxonomy_silva.make_otu_table_norm()
cluster.otu_cdhit.taxonomy_silva.make_plots()
cluster.otu_cdhit.taxonomy_silva.stats.nmds.run()
cluster.otu_cdhit.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_cdhit.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_cdhit.taxonomy_silva.comp_phyla.stats.nmds.run()
cluster.otu_cdhit.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_cdhit.taxonomy_silva.comp_tips.make_plots()
