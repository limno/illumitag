# Built-in modules #
from collections import defaultdict

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.clustering.composition import plots
from illumitag.common.csv_tables import CSVTable
from illumitag.common.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class Composition(object):
    """Base class for taxonomic compositioning"""

    all_paths = """
    /taxa_table.csv
    /graphs/
    """

    def __init__(self, parent, base_dir=None):
        # Parent #
        self.taxonomy, self.parent = parent, parent
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.composition_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__]
        # Taxa table #
        self.otu_csv = CSVTable(self.p.taxa_csv)

    @property_cached
    def frame(self):
        """Uses the otu_table and the taxonomy.assignments to make a new table"""
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.otu_table.iterrows():
            for otu_name, count in column.iteritems():
                taxa = self.assignments[otu_name]
                result[taxa[-1]][sample_name] += count
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Return #
        return result

    def make_taxa_table(self):
        """Convert to CSV"""
        self.taxa_table.to_csv(self.taxa_csv, sep='\t')