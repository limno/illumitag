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
from illumitag.fasta.single import FASTA
from illumitag.common.autopaths import DirectoryPath
from illumitag.clustering.taxonomy.crest import SimpleCrestTaxonomy
from illumitag.clustering.taxonomy.rdp import SimpleRdpTaxonomy

###############################################################################
# Get vars #
proj = illumitag.projects['evaluation']
pools = proj.pools
samples = [s for pool in proj for s in pool.samples]
cluster = illumitag.clustering.favorites.evaluation

# Load #
proj.load()

# Check bad samples #
s1 = illumitag.runs[1][3][31]
s2 = illumitag.runs[1][2][33]

# Compute Unifrac dist #
illumitag.helper.silvamod.amplify()
illumitag.helper.silvamod.amplified.graphs[1].plot(loglog=True)
cluster.otu_uparse.taxonomy_silva.make_filtered_centers()
cluster.otu_uparse.taxonomy_silva.stats.unifrac.run()
cluster.otu_uparse.taxonomy_silva.stats.unifrac.nmds.run()

# Make fraction graph #
proj.graphs[-1].plot()

# Get statistics #
proj.reporter.fraction_discarded

# Get clustering values #
r1, r2 = list(set([p.run for p in proj]))
r1.parse_report_xml()
r2.parse_report_xml()
print float(r1.report_stats['fwd']['DensityPF']) / float(r1.report_stats['fwd']['DensityRaw'])
print float(r2.report_stats['fwd']['DensityPF']) / float(r2.report_stats['fwd']['DensityRaw'])

# Check below 400 bp sequences #
folder = DirectoryPath(illumitag.projects['evaluation'].base_dir + "below_400/")
over = FASTA(folder + "reads.fasta")
def over_iterator(reads, max_length=400):
    for read in reads:
        if len(read) <= max_length: yield read
over.create()
for pool in pools: over.add_iterator(over_iterator(pool.good_barcodes.assembled.good_primers.qual_filtered))
over.close()
over.graphs[-1].plot()
crest = SimpleCrestTaxonomy(over, folder)
crest.assign()
crest.composition.graph.plot()
rdp = SimpleRdpTaxonomy(over, folder)
rdp.assign()
rdp.composition.graph.plot()
