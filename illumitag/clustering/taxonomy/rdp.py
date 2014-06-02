# Built-in modules #
import os

# Internal modules #
from illumitag.clustering.taxonomy import Taxonomy, SimpleTaxonomy
from illumitag.common.cache import property_cached
from illumitag.common.autopaths import AutoPaths
from illumitag.clustering.taxonomy import plots
from illumitag.fasta.single import FASTA
from illumitag.common.csv_tables import CSVTable
from illumitag.clustering.statistics import StatsOnTaxonomy

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #

###############################################################################
class RdpTaxonomy(Taxonomy):
    all_paths = """
    /otu_table.csv
    /graphs/
    /stats/
    /reads.taxonomy
    /composition.txt
    /assignment.txt
    """

    def __init__(self, fasta_path, parent, base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.rdp_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__]
        # Tables #
        self.otu_csv = CSVTable(self.p.otu_csv)
        # Composition tables #
        #self.comp_phyla = CompositionPhyla(self, self.p.comp_phyla)
        #self.comp_tips = CompositionTips(self, self.p.comp_tips)
        # Stats #
        self.stats = StatsOnTaxonomy(self)

    def assign(self):
        sh.rdp_multiclassifier('--conf=0.5','--hier_outfile='+self.p.composition, '--assign_outfile='+self.p.assignment)

    @property_cached
    def assignments(self):
        pass

###############################################################################
class SimpleRdpTaxonomy(SimpleTaxonomy):
    short_name = 'rdp'
    database = 'rdp'

    all_paths = """
    /composition.txt
    /assignments.txt
    """

    def assign(self):
        params = ['--hier_outfile='+self.p.composition, '--assign_outfile='+self.p.assignments]
        sh.rdp_multiclassifier(self.fasta, '--conf=0.5', *params)
        return self

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                line = line.strip('\n').split('\t')
                code = line.pop(0)
                species = tuple(x.strip('"') for x in line[1::3])
                result[code] = species
        return result