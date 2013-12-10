# Built-in modules #

# Internal modules #
import illumitag
from illumitag.common.cache import property_cached
from illumitag.common import natural_sort, prepend_to_file

# Third party modules #
import pandas

###############################################################################
class Taxonomy(object):
    """Can assign taxonomy to a FASTA file of 16S sequences."""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def at_level(self, level):
        return dict((k,v[level]) for k,v in self.assignments.items() if len(v) > level)

    def make_plots(self):
        for graph in self.graphs: graph.plot()

    @property_cached
    def otu_table(self):
        # Parse that custom output for the OTU table #
        result = pandas.DataFrame(self.otu.readmap.otu_sample_counts)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        # Remove unwanted #
        unwanted = ['Plastid', 'Mitochondrion', 'Thaumarchaeota', 'Crenarchaeota', 'Euryarchaeota']
        for otu_name in result:
            species = self.assignments[otu_name]
            if len(species) > 2 and species[2] in unwanted: result = result.drop(otu_name, 1)
        # Merge samples #
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
        # Return result #
        return result

    def make_otu_table(self):
        """Convert to CSV"""
        self.otu_table.to_csv(self.otu_csv, sep='\t')
        prepend_to_file(self.otu_csv, 'X')
