# Built-in modules #

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.analysis.otus import UclustOTUs, CdhitOTUs, OpenRefOTUs, StepOTUs
from illumitag.running.analysis_runner import AnalysisRunner
from illumitag.fasta.single import FASTA

# Third party modules #
from shell_command import shell_output

###############################################################################
class Analysis(object):
    """Analyzes an aggregate."""

    all_paths = """
    /reads/qiime_reads.fasta
    /otus/
    /logs/
    """

    def __repr__(self): return '<%s object on %s>' % \
                               (self.__class__.__name__, self.parent)
    def __iter__(self): return iter(self.pools)
    def __len__(self): return len(self.pools)
    def __getitem__(self, key): return self.pools[key]

    def __init__(self, aggregate):
        # Attributes #
        self.aggregate, self.parent = aggregate, aggregate
        # Inherited #
        self.pools, self.children = self.parent.pools, self.parent.pools
        self.base_dir = self.parent.base_dir
        self.meta_data_path = self.parent.meta_data_path
        # Dir #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner = AnalysisRunner(self)
        # FASTA #
        self.qiime_reads = FASTA(self.p.qiime_reads_fasta)
        # OTU picking #
        self.otu_uclust    = UclustOTUs(self)
        self.otu_cdhit     = CdhitOTUs(self)
        self.otu_open_ref  = OpenRefOTUs(self)
        self.otu_step      = StepOTUs(self)

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        self.runner.run_slurm()

    def combine_reads(self):
        paths = [pool.quality_reads.qiime_fasta.path for pool in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.qiime_reads.path))

    def run_uclust(self): self.otu_uclust.run()
    def run_cdhit(self): self.otu_cdhit.run()
    def run_open_ref(self): self.otu_open_ref.run()
    def run_progressive(self): self.otu_step.run()
