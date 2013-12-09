# Built-in modules #
import os

# Internal modules #
from illumitag.clustering.taxonomy import Taxonomy
from illumitag.common.cache import property_cached
from illumitag.common.autopaths import AutoPaths
from illumitag.clustering.taxonomy import plots
from illumitag.fasta.single import FASTA
from illumitag.common.csv_tables import CSVTable

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #

###############################################################################
class RdpTaxonomy(Taxonomy):
    all_paths = """
    /graphs/
    /taxa_table.csv
    /reads.taxonomy
    """

    def __init__(self, fasta_path, parent):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # Dir #
        self.base_dir = self.parent.p.rdp_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__]
        # Table #
        self.taxa_csv = CSVTable(self.p.taxa_csv)

    def assign(self):
        sh.rdp_multiclassifier('--conf=0.5','--hier_outfile='+self.composition_path,'--assign_outfile='+self.assignement_path)

    @property_cached
    def assignments(self):
        pass