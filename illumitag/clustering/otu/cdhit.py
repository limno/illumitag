# Built-in modules #
import os, shutil, re
from collections import defaultdict

# Internal modules #
import illumitag
from illumitag.common import natural_sort
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.common.cache import property_cached
from illumitag.common.tmpstuff import TmpFile
from illumitag.fasta.single import FASTA, FASTQ
from illumitag.clustering.otu import OTUs
from illumitag.clustering.taxonomy.crest import CrestTaxonomy

# Third party modules #
import sh, pandas
from shell_command import shell_output

# Constants #
home = os.environ['HOME'] + '/'
cdhit_script = home + 'share/cd-hit-otu-illumina/cd-hit-otu-all-single.pl'

###############################################################################
class CdhitOTUs(OTUs):
    """Will use cd-hit to create OTU clusters from a given FASTA file
    http://qiime.org/scripts/pick_otus.html"""

    short_name = 'cdhit'
    title = 'CD-HIT Illumina OTU picking'

    all_paths = """
    /all_reads.fastq
    /clusters/OTU.nr2nd.clstr
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
        # Main reads file here FASTQ #
        self.reads = FASTQ(self.p.all_reads)
        # Files #
        self.cdhit_clusters = FilePath(self.p.clstr)
        self.cdhit_centers = FASTA(self.p.clusters_dir + "OTU")
        self.centers = FASTA(self.p.centers)
        # Taxonomy #
        self.taxonomy_silva = CrestTaxonomy(self.centers, self, 'silvamod', self.p.silva)
        self.taxonomy_fw = CrestTaxonomy(self.centers, self, 'freshwater', self.p.fw_dir)
        # Preferred one #
        self.taxonomy = self.taxonomy_silva

    def run(self):
        # Combine reads but in fastq format this time #
        paths = [sample.renamed for sample in self.cluster]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))
        # Clean #
        shutil.rmtree(self.p.clusters_dir)
        # Run command #
        cdhit = sh.Command(cdhit_script)
        cdhit('-i', self.reads, '-o', self.p.clusters_dir, '-p', TmpFile.from_string('[ACTG]'))
        # Create the centers file with good names #
        self.cdhit_centers.rename_with_num('OTU_', self.centers)

    @property_cached
    def cluster_counts_table(self):
        """Create the unfiltered OTU table"""
        # Put results in a dict of dict #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for line in self.cdhit_clusters:
            if line.startswith('>'):
                otu = "OTU_%s" % line.split()[1]
                continue
            nums = re.findall(">run([0-9]+)_pool([0-9]+)_sample([0-9]+)_read([0-9]+)\.\.\.", line)
            if nums:
                run_num, pool_num, sample_num, read_num = map(int, nums[0])
                sample = illumitag.runs[run_num][pool_num-1][sample_num-1]
                name = sample.short_name
            else:
                nums = re.findall(">run([0-9]+)_sample([0-9]+)_read([0-9]+)\.\.\.", line)
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
        # Remove OTUs that are only one read #
        return result
