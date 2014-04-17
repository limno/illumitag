#!/usr/bin/env python2

"""
A script to run the clustering analyses on the evaluation samlpes.
"""

# Future #
from __future__ import division

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag

###############################################################################
# Get vars #
proj = illumitag.projects['evaluation']
samples = [s for pool in proj for s in pool.samples]
cluster = illumitag.clustering.favorites.evaluation

# Create samples #
for pool in proj: print pool