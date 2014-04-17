# Futures #
from __future__ import division

# Built-in modules #
from collections import Counter

# Internal modules #
import illumitag
from illumitag.common import natural_sort, prepend_to_file, flatten
from illumitag.common.cache import property_cached
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA

# Third party modules #
import pandas, numpy

###############################################################################
class Taxonomy(object):
    """Can assign taxonomy to a FASTA file of 16S sequences."""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def at_level(self, level):
        return dict((k,v[level]) for k,v in self.assignments.items() if len(v) > level)

    def make_plots(self):
        for graph in self.graphs: graph.plot()

    @property_cached
    def cluster_counts_table(self):
        """Parse that custom output for creating the unfiltered OTU table"""
        result = pandas.DataFrame(self.otu.readmap.otu_sample_counts)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        return result

    def make_otu_table(self):
        # Remove unwanted #
        result = self.cluster_counts_table.copy()
        unwanted = ['Plastid', 'Mitochondrion', 'Thaumarchaeota', 'Crenarchaeota', 'Euryarchaeota']
        for otu_name in result:
            species = self.assignments[otu_name]
            if len(species) > 2 and species[2] in unwanted: result = result.drop(otu_name, 1)
        # Merge samples when reruns #
        for name, data in result.iterrows():
            sample = [s for s in self.samples if s.short_name == name][0]
            rerun_row = result.loc[name]
            if not sample.info.get('rerun'): continue
            rn, p, n = sample.info['rerun']['run'], sample.info['rerun']['pool'], sample.info['rerun']['num']
            orig_sample = illumitag.runs[rn][p-1][n-1]
            orig_name = orig_sample.short_name
            orig_row = result.loc[orig_name]
            result.loc[orig_name] = orig_row + rerun_row
            result = result.drop(name)
        # Convert to CSV #
        result.to_csv(self.otu_csv, sep='\t')
        prepend_to_file(self.otu_csv, 'X')

    @property_cached
    def otu_table(self):
        return pandas.io.parsers.read_csv(self.otu_csv, sep='\t', index_col=0, encoding='utf-8')

    @property
    def otu_table_norm(self):
        """The same thing as otu_table but normalized so that the sum of a sample is always one"""
        return self.otu_table.apply(lambda x: x/x.sum(), axis=1).replace(numpy.inf, 0.0)

    def make_otu_table_norm(self):
        # Convert to CSV #
        self.otu_table_norm.to_csv(self.otu_csv_norm, sep='\t', float_format='%.5g')
        prepend_to_file(self.otu_csv_norm, 'X')

    def resample_otu_table(self, down_to=5000):
        # Eliminate samples that are under down_to #
        are_high = self.otu_table.sum(axis=1) > down_to
        old_frame = self.otu_table.loc[are_high,:]
        # Make a new table #
        new_frame = pandas.DataFrame()
        for i, sample in old_frame.iterrows():
            pop = flatten([[key]*val for key,val in freq.items()])
            smaller_pop = random.sample(pop, down_to)
            return Counter(smaller_pop)
            new_frame[sample] = sample

###############################################################################
class SimpleTaxonomy(object):

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.fasta)

    def __init__(self, fasta, base_dir):
        self.fasta = fasta if isinstance(fasta, FASTA) else FASTA(fasta)
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def phyla(self):
        return Counter(taxa[2] if len(taxa) > 2 else taxa[-1] for taxa in self.assignments.values())