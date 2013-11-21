# Built-in modules #

# Internal modules #
import illumitag
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.running.cluster_runner import ClusterRunner
from illumitag.clustering.uparse import UparseOTUs

# Third party modules #
from shell_command import shell_output

###############################################################################
class Cluster(object):
    """Analyzes a group of samples."""

    all_paths = """
    /reads/all_reads.fasta
    /otus/
    /logs/
    """

    def __repr__(self): return '<%s object with %i samples>' % (self.__class__.__name__, len(self.samples))
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key): return self.samples[key]

    def __init__(self, samples, name):
        # Save samples #
        self.name = name
        self.samples, self.children = samples, samples
        # Figure out pools #
        self.pools = list(set([s.pool for s in self.samples]))
        # Load them #
        for p in self.pools: p.load()
        # Dir #
        self.base_dir = illumitag.view_dir + "clusters/" + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner = ClusterRunner(self)
        # FASTA #
        self.reads = FASTA(self.p.all_reads_fasta)
        # OTU picking #
        self.otu_uparse = UparseOTUs(self)

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        self.runner.run_slurm()

    def combine_reads(self):
        paths = [sample.fasta.path for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))

    def run_uparse(self): self.otu_uparse.run()