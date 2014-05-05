# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.fasta.single import FASTA
from illumitag.common.autopaths import AutoPaths
from illumitag.common.slurm import nr_threads
from illumitag.common.cache import property_cached
from illumitag.common.csv_tables import CSVTable
from illumitag.clustering.statistics import StatsOnTaxonomy
from illumitag.clustering.taxonomy import Taxonomy, SimpleTaxonomy
from illumitag.clustering.taxonomy import plots
from illumitag.clustering.composition.phyla import CompositionPhyla
from illumitag.clustering.composition.tips import CompositionTips

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #
databases = {
    'silvamod': home + "share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta",
    'freshwater': home + "share/LCAClassifier/parts/flatdb/freshwater/freshwater.fasta"
}

###############################################################################
class CrestTaxonomy(Taxonomy):
    all_paths = """
    /otu_table.csv
    /otu_table_norm.csv
    /centers.fasta
    /graphs/
    /stats/
    /comp_phyla/
    /comp_tips/
    /crest_hits.xml
    /crest_composition.txt
    /crest_tree.txt
    /crest_assignments.txt
    """

    def __init__(self, fasta_path, parent, database='silvamod', base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # The database to use #
        self.database = database
        self.database_path = databases[database]
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.crest_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__]
        # OTU table #
        self.otu_csv = CSVTable(self.p.otu_csv)
        self.otu_csv_norm = CSVTable(self.p.otu_csv_norm)
        # Filtered centers file #
        self.centers = FASTA(self.p.centers)
        # Composition tables #
        self.comp_phyla = CompositionPhyla(self, self.p.comp_phyla)
        self.comp_tips = CompositionTips(self, self.p.comp_tips)
        # Stats #
        self.stats = StatsOnTaxonomy(self)

    def assign(self):
        # Run #
        sh.megablast('-a', nr_threads, '-i', self.fasta, '-d', self.database_path, '-b100', '-v100', '-m7', '-o', self.p.crest_hits)
        if os.path.getsize(self.p.crest_hits) == 0: raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        sh.classify(self.p.crest_hits, '-p', '-o', '-d', self.database)
        shutil.move(self.p.crest_hits[:-4] + '_Composition.txt', self.p.crest_composition)
        shutil.move(self.p.crest_hits[:-4] + '_Tree.txt', self.p.crest_tree)
        shutil.move(self.p.crest_hits[:-4] + '_Assignments.txt', self.p.crest_assignments)
        # Clean up #
        comp_path = home + 'ILLUMITAG/All_Composition.txt'
        if os.path.exists(comp_path): os.remove(comp_path)
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")
        # Return #
        return self

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                code, species = line.split('\t')
                result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result

###############################################################################
class SimpleCrestTaxonomy(SimpleTaxonomy, CrestTaxonomy):
    short_name = 'crest'

    def __init__(self, fasta, base_dir, database='silvamod'):
        SimpleTaxonomy.__init__(self, fasta, base_dir)
        self.database = database
        self.database_path = databases[database]
