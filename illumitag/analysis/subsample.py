# Built-in modules #
import random
from collections import Counter

# Internal modules #
from illumitag.analysis.otus import OTUs
from illumitag.common.autopaths import AutoPaths
from illumitag.common.csv_tables import TSVTable
from illumitag.analysis.statistics import StatsOnOTU
from illumitag.graphs import otu_plots

# Third party modules #
import pandas

# Constants #

###############################################################################
class SubsampledOTUs(OTUs):
    dist_method = 'bray'

    #all_paths = OTUs.all_paths + """
    #/table/subsampled.csv
    #"""

    def __repr__(self): return '<%s object of %s>' % \
                               (self.__class__.__name__, self.parent)

    def __init__(self, parent):
        # Save parent #
        self.full_otu, self.parent = parent, parent
        # Names #
        self.short_name = parent.short_name + '_subsampled'
        self.method = parent.method + ' (Subsampled)'
        # Inherited #
        self.pools = self.parent.pools
        self.qiime_reads = self.parent.qiime_reads
        self.meta_data_path = self.parent.meta_data_path
        # Paths #
        self.base_dir = self.parent.p.subsampled_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Other #
        self.taxonomy = None
        # Files #
        self.table = TSVTable(self.p.csv_table)
        self.table_filtered = TSVTable(self.p.csv_table_filtered)
        self.table_transposed = TSVTable(self.p.csv_table_transposed)
        # Children #
        self.stats = StatsOnOTU(self, self.table_filtered)

    def run(self, down_to=None):
        self.subsample()
        self.make_otu_plots()
        self.compute_stats()

    def subsample(self, down_to=None):
        # Parse #
        otus = pandas.read_csv(self.full_otu.table_filtered.path, sep = '\t', index_col=0)
        # Determine down_to #
        sums = otus.sum(axis=1)
        if not down_to: self.down_to = min(sums)
        else:
            self.down_to = down_to
            otus = otus.drop(sums[sums < self.down_to].keys())
        # Empty frame #
        subotus = pandas.DataFrame(columns=otus.columns, index=otus.index, dtype=int)
        # Do it #
        for sample_name in otus.index:
            row = otus.loc[sample_name]
            weighted_choices = list(row[row != 0].iteritems())
            population = [val for val, count in weighted_choices for i in range(count)]
            sub_pop = random.sample(population, self.down_to)
            frequencies = Counter(sub_pop)
            new_row = pandas.Series(frequencies.values(), index=frequencies.keys(), dtype=int)
            subotus.loc[sample_name] = new_row
        # Output it #
        subotus.to_csv(self.table.path, sep='\t', na_rep='0')
        self.table.to_integer(path=self.table_filtered.path)
        # Transposed #
        self.table_filtered.transpose(path=self.table.path)
        self.table_filtered.transpose(path=self.table_transposed.path)

    def make_otu_plots(self):
        for cls_name in otu_plots.__all__[1:]:
            cls = getattr(otu_plots, cls_name)
            cls(self).plot()
