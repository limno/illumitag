# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.cluster.taxonomy import Taxonomy
from illumitag.common.slurm import nr_threads
from illumitag.common.cache import property_cached

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #
train_db_path = home + 'glob/16s/trainset/trainset9_032012.pds.fasta'
train_tax_path = home + 'glob/16s/trainset/trainset9_032012.pds.tax'
silva_db_path = home + "/share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta"

###############################################################################
class CrestTaxonomy(Taxonomy):
    all_paths = """
    /reads.taxonomy
    /silva_hits.xml
    /silva_composition.txt
    /silva_tree.txt
    /silva_assignments.txt
    """

    def assign(self):
        # Run #
        sh.megablast('-a', nr_threads, '-i', self.fasta, '-d', silva_db_path, '-b100', '-v100', '-m7', '-o', self.p.silva_hits)
        if os.path.getsize(self.p.silva_hits) == 0: raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        sh.classify(self.p.silva_hits, '-p', '-o', '-d', 'silvamod')
        shutil.move(self.p.silva_hits[:-4] + '_Composition.txt', self.p.silva_composition)
        shutil.move(self.p.silva_hits[:-4] + '_Tree.txt', self.p.silva_tree)
        shutil.move(self.p.silva_hits[:-4] + '_Assignments.txt', self.p.silva_assignments)

    @property_cached
    def assignments(self):
        result = {}
        for line in open(self.p.assignments):
            code, species = line.split('\t')
            result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result
