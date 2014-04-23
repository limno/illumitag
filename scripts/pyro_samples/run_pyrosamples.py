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
for sample in illumitag.pyrosamples: sample.report_loss()
for sample in illumitag.pyrosamples: sample.process()
for sample in illumitag.pyrosamples: sample.raw_fastq.fastqc()
for sample in illumitag.pyrosamples: sample.fasta.graphs[-1].plot()

# Cluster #
cluster = illumitag.clustering.favorites.pyro_comparison
cluster.combine_reads()
cluster.set_size(400)

# UParse #
cluster.otu_uparse.taxonomy_silva.comp_phyla.graphs[-1].plot()

# UClust #
pass

# CD-hit #
for sample in illumitag.pyrosamples: sample.make_fastq()

# A graph #
cluster.reads.graphs[-1].plot()