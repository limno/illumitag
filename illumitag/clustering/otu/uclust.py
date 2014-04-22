# Built-in modules #
import shutil, re
from collections import defaultdict

# Internal modules #
import illumitag
from illumitag.common import natural_sort
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.common.cache import property_cached
from illumitag.fasta.single import FASTA
from illumitag.clustering.otu import OTUs
from illumitag.clustering.taxonomy.crest import CrestTaxonomy

# Third party modules #
import sh, pandas

###############################################################################
class UclustOTUs(OTUs):
    """Will use uclust via the qimme wraper to create OTU clusters from a given FASTA file
    http://qiime.org/scripts/pick_otus.html"""

    short_name = 'uclust'
    title = 'UCLUST-QIIME denovo picking'

    all_paths = """
    /clusters/clusters.uc
    /clusters/qiime.log
    /clusters/all_otus.txt
    /clusters/all_centers.fasta
    /centers.fasta
    /otus.txt
    /taxonomy_silva/
    /taxonomy_fw/
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
        self.all_otus = FilePath(self.p.all_otus)
        self.all_centers = FASTA(self.p.all_centers)
        self.otus = FilePath(self.base_dir + "otus.txt")
        self.centers = FASTA(self.base_dir + "centers.fasta")
        # Taxonomy #
        self.taxonomy_silva = CrestTaxonomy(self.centers, self, 'silvamod', self.p.silva)
        self.taxonomy_fw = CrestTaxonomy(self.centers, self, 'freshwater', self.p.fw_dir)
        # Preferred one #
        self.taxonomy = self.taxonomy_silva

    def run(self):
        # Clean #
        shutil.rmtree(self.p.clusters_dir)
        # Run command #
        pick_otus = sh.Command('pick_otus.py')
        pick_otus('-m', 'uclust', '-i', self.reads, '-o', self.p.clusters_dir)
        # Move into place #
        base_name = self.p.clusters_dir + self.reads.prefix
        shutil.move(base_name + '_otus.txt', self.all_otus)
        shutil.move(base_name + '_otus.log', self.p.qiime_log)
        shutil.move(base_name + '_clusters.uc', self.p.clusters_uc)
        # Remove OTUs that are only one read #
        def filter_singletons(f):
            for line in f:
                line = line.split()
                if len(line) > 2: yield '\t'.join(line) + '\n'
        self.otus.writelines(filter_singletons(self.all_otus))
        # Create the centers file that is missing #
        pick_rep = sh.Command('pick_rep_set.py')
        pick_rep('-i', self.all_otus, '-f', self.reads, '-o', self.all_centers)
        # Remake the centers file without the filtered OTUs #
        self.otus_to_keep = [line.split()[0] for line in self.otus]
        def filter_otus(f):
            for seq in f:
                if seq.id in self.otus_to_keep: yield seq
        self.centers.write(filter_otus(self.all_centers))

    @property_cached
    def cluster_counts_table(self):
        """Create the unfiltered OTU table"""
        # Put results in a dict of dicts #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for line in self.otus:
            # Parse the line #
            contents = line.split()
            otu, reads = contents[0], contents[1:]
            # Parse the hits #
            for r in reads:
                nums = re.findall("run([0-9]+)_pool([0-9]+)_sample([0-9]+)_read([0-9]+)", r)
                if nums:
                    run_num, pool_num, sample_num, read_num = map(int, nums[0])
                    sample = illumitag.runs[run_num][pool_num-1][sample_num-1]
                    name = sample.short_name
                else:
                    nums = re.findall("run([0-9]+)_sample([0-9]+)_read([0-9]+)", r)
                    run_num, sample_num, read_num = map(int, nums[0])
                    sample = [s for s in illumitag.presamples+illumitag.pyrosamples if s.run_num==run_num and s.num==sample_num][0]
                    name = sample.short_name
                # Count #
                result[otu][name] += 1
        # Return #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis=1)
        return result