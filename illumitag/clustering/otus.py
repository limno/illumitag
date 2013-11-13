# Built-in modules #
import os, re
from collections import defaultdict

# Internal modules #
import illumitag
from illumitag.common import natural_sort
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.common.csv_tables import CSVTable
from illumitag.fasta.single import FASTA, SizesFASTA
from illumitag.clustering import otu_plots

# Third party modules #
import sh, pandas

# Constants #
home = os.environ['HOME'] + '/'
green_genes_db_path = home + "share/green_genes/rep_set/97_otus.fasta"

###############################################################################
class UparseOTUs(object):
    """Will use uparse to create an OTU table from a given FASTA file
    http://www.nature.com/doifinder/10.1038/nmeth.2604"""

    short_name = 'uparse'

    all_paths = """
    /derep.fasta
    /sorted.fasta
    /centers.fasta
    /readmap.uc
    /table.csv
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
        self.derep = SizesFASTA(self.p.derep)
        self.sorted = SizesFASTA(self.p.sorted)
        self.centers = FASTA(self.p.centers)
        self.readmap = UClusterFile(self.p.readmap)
        self.table = CSVTable(self.p.table)
        # Graphs #
        self.graphs = [getattr(otu_plots, cls_name)(self) for cls_name in otu_plots.__all__]

    def checks(self):
        assert len(self.reads) == len(self.derep)
        assert len(self.reads) == len(self.readmap)

    def run(self):
        ## Dereplicate #
        #sh.usearch7("--derep_fulllength", self.reads, '-output', self.derep, '-sizeout')
        ## Order by size and kill singeltons #
        #sh.usearch7("--sortbysize", self.derep, '-output', self.sorted, '-minsize', 2)
        ## Compute the centers #
        #sh.usearch7("--cluster_otus", self.sorted, '-otus', self.centers)
        ## Rename the centers #
        #self.centers.rename_with_num('OTU_')
        ## Map the reads back to the centers #
        #sh.usearch7("-usearch_global", self.reads, '-db', self.centers, '-strand', 'plus', '-id', 0.97, '-uc', self.readmap)
        # Parse that custom output for the OTU table #
        self.frame = pandas.DataFrame(self.readmap.otu_sample_counts)
        self.frame = self.frame.fillna(0)
        self.frame = self.frame.reindex_axis(sorted(self.frame.columns, key=natural_sort), axis=1)
        # Convert to CSV #
        self.frame.to_csv(self.table, sep='\t')

###############################################################################
class UClusterFile(FilePath):
    """A special format outputed by uparse
    An example line:
    H       1474    422     97.6    +       0       0       422M    run4_pool1_sample1_read2        OTU_1474

    The corresponding legend:
    # 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
    """

    def __len__(self): return self.count
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0

    @property
    def count(self):
        return int(sh.wc('-l', self.path).stdout.split()[0])

    @property
    def otu_sample_counts(self):
        # Put results in a dict of dict #
        result = defaultdict(lambda: defaultdict(int))
        # Loop #
        for line in self:
            # Parse the line #
            kind, num, length, similarity, strand, start, seed, alignment, query, target = line.split()
            # Skip no hits #
            if target == '*': continue
            # Parse the hit #
            nums = re.findall("run([0-9]+)_pool([0-9]+)_sample([0-9]+)_read([0-9]+)", query)
            run_num, pool_num, sample_num, read_num = map(int, nums[0])
            # Get the object #
            sample = illumitag.runs[run_num][pool_num-1][sample_num-1]
            # Count #
            result[target][sample.short_name] += 1
        # Return #
        return result