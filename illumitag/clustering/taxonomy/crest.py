# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.clustering.taxonomy import Taxonomy
from illumitag.clustering.taxonomy import plots
from illumitag.common.slurm import nr_threads
from illumitag.common.cache import property_cached
from illumitag.common.csv_tables import CSVTable

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #
silva_db_path = home + "/share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta"
fw_db_path = home + "/share/LCAClassifier/parts/flatdb/freshwater/freshwater.fasta"

###############################################################################
class CrestTaxonomy(Taxonomy):
    all_paths = """
    /graphs/
    /taxa_table.csv
    /reads.taxonomy
    /silva_hits.xml
    /silva_composition.txt
    /silva_tree.txt
    /silva_assignments.txt
    """

    def __init__(self, fasta_path, parent, database=None, base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # The database to use #
        if database is None:      self.database = silva_db_path
        elif database == 'silva': self.database = silva_db_path
        elif database == 'fw':    self.database = fw_db_path
        else:                     self.database = database
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.crest_dir
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graphs #
        self.graphs = [getattr(plots, cls_name)(self) for cls_name in plots.__all__]
        # Table #
        self.taxa_csv = CSVTable(self.p.taxa_csv)

    def assign(self):
        # Run #
        sh.megablast('-a', nr_threads, '-i', self.fasta, '-d', self.database, '-b100', '-v100', '-m7', '-o', self.p.silva_hits)
        if os.path.getsize(self.p.silva_hits) == 0: raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        sh.classify(self.p.silva_hits, '-p', '-o', '-d', 'silvamod')
        shutil.move(self.p.silva_hits[:-4] + '_Composition.txt', self.p.silva_composition)
        shutil.move(self.p.silva_hits[:-4] + '_Tree.txt', self.p.silva_tree)
        shutil.move(self.p.silva_hits[:-4] + '_Assignments.txt', self.p.silva_assignments)

    @property_cached
    def assignments(self):
        result = {}
        with open(self.p.assignments, 'r') as handle:
            for line in handle:
                code, species = line.split('\t')
                result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result
