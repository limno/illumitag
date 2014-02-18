#!/usr/bin/env python2

"""
A script to run the pyrosamples.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import illumitag

###############################################################################
illumitag.demultiplexer.run()
illumitag.pyrosamples[0].extract()
illumitag.pyrosamples[0].clean()
illumitag.pyrosamples[0].raw_fastq.fastqc()

execfile("/home/lucass/repos/illumitag/scripts/single_use/pyro_plots.py")
cluster.otu_uparse.taxonomy_silva.comp_phyla.graphs[-1].plot()

# Check size distribution #
for s in illumitag.pyrosamples: s.