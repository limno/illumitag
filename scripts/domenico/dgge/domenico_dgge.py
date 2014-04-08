#!/usr/bin/env python2
#-*- coding: utf-8 -*-

"""
A custom made script to run some sequence similarity searches on some OTUs.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ipython -i domenico_dgge.py records.fasta
"""

###############################################################################
# The libraries we need #
import inspect, sys, os
import illumitag, parallelblast
# Check usage #
if len(sys.argv) != 2: sys.exit(sys.modules[__name__].__doc__)
# Get the shell arguments #
records_path = os.path.expanduser(sys.argv[1])
# Check that the path is valid #
if not os.path.exists(records_path): raise Exception("No file at %s." % records_path)

# Get this directory #
file_name = inspect.getframeinfo(inspect.currentframe()).filename
this_dir = os.path.dirname(os.path.abspath(file_name)) + '/'
# Copy OTU sequences here #
centers_path = this_dir + 'centers.fasta'
cluster = illumitag.clustering.favorites.danube
cluster.otu_uparse.centers.copy(centers_path)
# Make database #
db = parallelblast.BLASTdb(centers_path)
db.makeblastdb()

# Search #
params = {'executable': "~/share/blastplus/blastn",
          '-outfmt': 0,
          '-evalue': 1e-2,
          '-perc_identity': 95,
          '-num_alignments': 0,
          '-num_threads': 16}
search = parallelblast.BLASTquery(records_path, db, params)
search.run()