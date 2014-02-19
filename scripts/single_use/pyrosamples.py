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
for sample in illumitag.pyrosamples: sample.extract()
for sample in illumitag.pyrosamples: sample.clean()
for sample in illumitag.pyrosamples: sample.process()
for sample in illumitag.pyrosamples: sample.raw_fastq.fastqc()
for sample in illumitag.pyrosamples: sample.fasta.graphs[-1].plot()

execfile("/home/lucass/repos/illumitag/scripts/single_use/pyro_plots.py")
cluster.otu_uparse.taxonomy_silva.comp_phyla.graphs[-1].plot()
cluster.set_size(400)
cluster.reads.graphs[-1].plot()
