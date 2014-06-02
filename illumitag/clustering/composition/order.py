# Built-in modules #
from collections import defaultdict

# Internal modules #
from illumitag.clustering.composition import Composition
from illumitag.common.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class CompositionOrder(Composition):
    """The taxa table contains only the order level, and sometimes
    the class levels, when phyla are over a certain threshold.
    Additionally, low abundance taxa are grouped into 'Others'"""

    @property_cached
    def taxa_table(self):
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.parent.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.parent.assignments[otu_name]
                order = taxa[4] if len(taxa) > 4 else taxa[-1]
                result[order][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Ungroup high abundance #
        high_abundance = result.sum() > 300000
        for order in [k for k,v in high_abundance.iteritems() if v]:
            # Find OTUs #
            cond = lambda taxa: len(taxa) > 4 and taxa[4] == order
            otus = [otu for otu,taxa in self.parent.assignments.items() if cond(taxa)]
            # Make new columns #
            new_columns = defaultdict(lambda: defaultdict(int))
            for sample_name, column in self.parent.otu_table.iterrows():
                for otu_name in otus:
                    count = column[otu_name]
                    taxa = self.parent.assignments[otu_name]
                    family = taxa[5] if len(taxa) > 5 else taxa[-1]
                    if not family.endswith('(family)'): family += ' (family)'
                    new_columns[family][sample_name] += count
            new_columns = pandas.DataFrame(new_columns)
            # Check #
            assert (new_columns.sum(axis=1) == result[order]).all()
            # Switch them in place #
            result = result.drop(order, axis=1)
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
