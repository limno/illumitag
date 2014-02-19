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

# Create the cluster #
cluster = illumitag.clustering.Cluster(illumitag.runs[0][0].samples.children, 'test')

# Run it #
for s in cluster: s.process()
cluster.combine_reads()
cluster.set_size()
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()
cluster.otu_uparse.taxonomy_silva.make_otu_table()
cluster.otu_uparse.taxonomy_silva.make_plots()
cluster.otu_uparse.taxonomy_silva.stats.nmds.run()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_taxa_table()
cluster.otu_uparse.taxonomy_silva.comp_tips.make_plots()

###############################################################################
# Pyro cluster #
samples = [s for s in illumitag.pyrosamples]
samples += illumitag.runs[1][0][0:8]
samples += illumitag.runs[1][1][0:8]
samples += illumitag.runs[2][0][0:8]
samples += illumitag.runs[3][7][30:33]
samples += illumitag.runs[3][7][13:16]
samples += illumitag.runs[3][7][43:46]
cluster = illumitag.clustering.Cluster(samples, 'pyro_comparison')

# Domenico #
samples = [s for pool in illumitag.runs[4][0:3] for s in pool.samples if s.used]
samples += illumitag.runs[4][3][0:13]
cluster = illumitag.clustering.Cluster(samples, 'domenico')

# Mixed eval #
samples = [s for s in illumitag.presamples]
samples += illumitag.runs[1][0][0:8]
samples += illumitag.runs[1][1][0:8]
samples += illumitag.runs[2][0][0:8]
cluster = illumitag.clustering.Cluster(samples, 'mixed_evaluation')

# Other clusters #
cluster = illumitag.clustering.Cluster(illumitag.presamples, 'new_lab_test_with')

# Soda lakes #
temporal = [s for s in illumitag.runs[3][7].samples if s.used]
spatial = [s for s in illumitag.runs[4][3].samples if s.used and s.group_name != 'JDS_2007']
samples = temporal + spatial
cluster = illumitag.clustering.Cluster(samples, 'soda')

###############################################################################
# Inga's cluster #
samples = [s for s in illumitag.runs[3][6].samples if s.used]
cluster = illumitag.clustering.Cluster(samples, 'inga')
for s in tqdm(cluster): s.process()
cluster.combine_reads()

# Anna's cluster #
from tqdm import tqdm
samples = [s for s in illumitag.runs[4][4].samples if s.used]
samples += [s for s in illumitag.runs[4][5].samples if s.used]
cluster = illumitag.clustering.Cluster(samples, 'anna')
for s in tqdm(cluster): s.process()
cluster.combine_reads()
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()

# Jerome's cluster #
from tqdm import tqdm
samples = [s for s in illumitag.runs[4][6].samples if s.used]
samples += [s for s in illumitag.runs[4][7].samples if s.used]
cluster = illumitag.clustering.Cluster(samples, 'jerome')
for s in tqdm(cluster): s.process()
cluster.combine_reads()
cluster.otu_uparse.run()
cluster.otu_uparse.taxonomy_silva.assign()

# Monica's cluster #
from tqdm import tqdm
illumitag.runs[5][3].create_samples()
illumitag.runs[5][4].create_samples()
illumitag.runs[5][5].create_samples()
samples =  [s for s in illumitag.runs[5][3].samples if s.used]
samples += [s for s in illumitag.runs[5][4].samples if s.used]
samples += [s for s in illumitag.runs[5][5][0:11] if s.used]
cluster = illumitag.clustering.Cluster(samples, 'monica')
for s in tqdm(cluster): s.process()
