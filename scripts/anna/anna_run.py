#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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
# Get the cluster #
cluster = illumitag.clustering.favorites.anna

# Stuff #
cluster.otu_uparse.taxonomy_silva.make_filtered_centers()

# Run seqenv #
cluster.run(steps=[{'otu_uparse.seqenv.run': dict(threads=False)}])
cluster.run_slurm(steps=[{'otu_uparse.seqenv.run': dict(threads=False)}])