# Built-in modules #
from collections import defaultdict

# Internal modules #
from illumitag.clustering.composition import Composition
from illumitag.common.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class CompositionPhyla(Composition):
    """Base class for taxonomical compositioning"""

    @property_cached
    def taxa_table(self):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.parent.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.parent.assignments[otu_name]
                phyla = taxa[2] if len(taxa) > 2 else taxa[-1]
                result[phyla][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Ungroup high abundance #
        high_abundance = result.sum() > 300000
        for phyla in [k for k,v in high_abundance.iteritems() if v]:
            # Find OTUs #
            cond = lambda taxa: len(taxa) > 2 and taxa[2] == phyla
            otus = [otu for otu,taxa in self.parent.assignments.items() if cond(taxa)]
            # Make new columns #
            new_columns = defaultdict(lambda: defaultdict(int))
            for sample_name, column in self.parent.otu_table.iterrows():
                for otu_name in otus:
                    count = column[otu_name]
                    taxa = self.parent.assignments[otu_name]
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
