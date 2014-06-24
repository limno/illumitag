#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run the clustering analyses for Silke's experiment.
"""

# Future #
from __future__ import division

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag

###############################################################################
# Get the cluster #
cluster = illumitag.clustering.favorites.anna

# Run UPARSE with different threshold #
cluster.otu_uparse.run(threshold=1.0)

# Other stuff #
cluster.otu_uparse.taxonomy_silva.assign()
cluster.otu_uparse.taxonomy_silva.make_otu_table()
cluster.otu_uparse.taxonomy_silva.make_otu_table_norm()
cluster.otu_uparse.taxonomy_silva.make_plots()
cluster.otu_uparse.taxonomy_silva.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.make_filtered_centers()

# Run seqenv #
cluster.otu_uparse.seqenv.run(threshold=1.0)

# Run seqenv via SLURM #
cluster.run(steps=[{'otu_uparse.seqenv.run': dict(threads=False)}])
cluster.run_slurm(steps=[{'otu_uparse.seqenv.run': dict(threads=False)}], time="1-00:00:00")
"bash -x /home/lucass/share/seqenv/SEQenv_v0.8/SEQenv_samples.sh -f /home/lucass/ILLUMITAG/views/clusters/anna/otus/uparse/taxonomy_silva/centers.fasta -s /home/lucass/ILLUMITAG/views/clusters/anna/otus/uparse/seqenv/abundances.csv -n 1000 -p -c 16"