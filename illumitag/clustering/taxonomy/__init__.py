# Built-in modules #
from collections import defaultdict

# Internal modules #
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
    def taxa_table(self):
        """Uses the otu_table and the taxonomy.assignments to make a new table"""
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.assignments[otu_name]
                phyla = taxa[2] if len(taxa) > 2 else taxa[-1]
                result[phyla][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Condition #
        if len(result.columns) > 12: pass
        # Ungroup high abundance #
        high_abundance = result.sum() > 300000
        for phyla in [k for k,v in high_abundance.iteritems() if v]:
            # Find OTUs #
            cond = lambda taxa: len(taxa) > 2 and taxa[2] == phyla
            otus = [otu for otu,taxa in self.assignments.items() if cond(taxa)]
            # Make new columns #
            new_columns = defaultdict(lambda: defaultdict(int))
            for sample_name, column in self.otu_table.iterrows():
                for otu_name in otus:
                    count = column[otu_name]
                    taxa = self.assignments[otu_name]
                    clss = taxa[3] if len(taxa) > 3 else taxa[-1]
                    if not clss.endswith('(class)'): clss += ' (class)'
                    new_columns[clss][sample_name] += count
            new_columns = pandas.DataFrame(new_columns)
            # Check #
            assert (new_columns.sum(axis=1) == result[phyla]).all()
            # Switch them in place #
            result = result.drop(phyla, axis=1)
            result = result.join(new_columns)
        # Group low abundant into 'others' #
        low_abundance = result.sum() < 3000
        other_count = result.loc[:, low_abundance].sum(axis=1)
        result = result.loc[:, ~low_abundance]
        result['Others'] = other_count
        # Sort the table by sum #
        sums = result.sum()
        sums.sort(ascending=False)
        result = result.reindex_axis(sums.keys(), axis=1)
        # Return result #
        return result

    def make_taxa_table(self):
        """Convert to CSV"""
        self.taxa_table.to_csv(self.taxa_csv, sep='\t')

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
        # Return result #
        return result

    def make_otu_table(self):
        """Convert to CSV"""
        self.otu_table.to_csv(self.otu_csv, sep='\t')
        prepend_to_file(self.otu_csv, 'X')
