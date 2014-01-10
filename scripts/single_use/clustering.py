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
# Create the samples #
illumitag.runs[4][0].create_samples()

# Make a cluster of samples #
samples = [s for pool in illumitag.runs[4][0:3] for s in pool.samples if s.used]
samples += illumitag.runs[4][3][0:13]
cluster = illumitag.clustering.Cluster(samples, 'domenico')
for s in cluster: s.process()

# Mixed eval #
samples = [s for s in illumitag.presamples]
samples += illumitag.runs[1][0][0:8]
samples += illumitag.runs[1][1][0:8]
samples += illumitag.runs[2][0][0:8]
cluster = illumitag.clustering.Cluster(samples, 'mixed_evaluation')

# Other clusters #
cluster = illumitag.clustering.Cluster(illumitag.runs[0][0].samples.children, 'test')
cluster = illumitag.clustering.Cluster(illumitag.presamples, 'new_lab_test_with')

# Run it #
cluster.combine_reads()
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()
cluster.otu_uparse.taxonomy_silva.make_otu_table()
cluster.otu_uparse.taxonomy_silva.make_plots()
cluster.otu_uparse.taxonomy_silva.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_plots()

# Soda lakes #
illumitag.runs[3][7].run_slurm()
illumitag.runs[4][3].run_slurm()
temporal = [s for s in illumitag.runs[3][7].samples if s.used]
spatial = [s for s in illumitag.runs[4][3].samples if s.used and s.group_name != 'JDS_2007']
samples = temporal + spatial
cluster = illumitag.clustering.Cluster(samples, 'soda')
for s in cluster: s.process()
cluster.combine_reads()
cluster.export_metadata()

###############################################################################
# Inga's cluster #
illumitag.runs[3][6].create_samples()
samples = [s for s in illumitag.runs[3][6].samples if s.used]
cluster = illumitag.clustering.Cluster(samples, 'inga')
for s in tqdm(cluster): s.process()
cluster.combine_reads()
