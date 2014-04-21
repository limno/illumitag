# Built-in modules #
import shutil

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA, SizesFASTA
from illumitag.clustering.otu import OTUs
from illumitag.clustering.taxonomy.crest import CrestTaxonomy

# Third party modules #
import sh

###############################################################################
class UclustOTUs(OTUs):
    """Will use uclust via the qimme wraper to create OTU clusters from a given FASTA file
    http://qiime.org/scripts/pick_otus.html"""

    short_name = 'uclust'
    title = 'UCLUST-QIIME denovo picking'

    all_paths = """
    /uclust_picked_otus/otus.txt
    /uclust_picked_otus/otus.log
    /uclust_picked_otus/clusters.uc
    /taxonomy_silva/
    /graphs/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, cluster):
        # Save parent #
        self.cluster, self.parent = cluster, cluster
        # Inherited #
        self.samples = self.parent.samples
        # Paths #
        self.base_dir = self.parent.p.otus_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Main FASTA file #
        self.reads = self.parent.reads
        # Files #
        #self.derep = SizesFASTA(self.p.derep)
        #self.sorted = SizesFASTA(self.p.sorted)
        #self.centers = FASTA(self.p.centers)
        # Taxonomy #
        #self.taxonomy_silva = CrestTaxonomy(self.centers, self, 'silvamod', self.p.silva)
        # Preferred one #
        #self.taxonomy = self.taxonomy_silva

    def run(self):
        # Clean #
        shutil.rmtree(self.p.picked_dir)
        # Prepare #
        pick_otus = sh.Command('pick_otus.py')
        # Run command #
        pick_otus('-m', 'uclust', '-i', self.reads)
        # Move into place #
        base_name = self.p.clusters_dir + self.qiime_reads.prefix
        shutil.move(base_name + '_otus.txt', self.p.clusters_otus_txt)
        shutil.move(base_name + '_otus.log', self.p.clusters_otus_log)
        shutil.move(base_name + '_clusters.uc', self.p.clusters_uc)
