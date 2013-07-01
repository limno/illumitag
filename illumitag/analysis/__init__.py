# Built-in modules #

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.analysis.otus import DenovoOTUs, OpenRefOTUs, StepOTUs
from illumitag.analysis.subsample import SubsampledOTUs
from illumitag.running.analysis_runner import AnalysisRunner
from illumitag.common.slurm import SLURMJob
from illumitag.fasta.single import FASTA

# Third party modules #
from shell_command import shell_output

###############################################################################
class Analysis(object):
    """Analyzes an aggregate."""

    all_paths = """
    /reads/qiime_reads.fasta
    /otus/
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
        self.otu_denovo    = DenovoOTUs(self)
        self.otu_open_ref  = OpenRefOTUs(self)
        self.otu_step      = StepOTUs(self)
        self.otu_subsample = SubsampledOTUs(self)

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, steps=None, **kwargs):
        # Check loaded #
        if not self.loaded: self.load()
        # Make script #
        command = """steps = %s
                     analysis = [pj for p in illumitag.pools if str(pj)=='%s'][0]
                     analysis(steps)""" % (steps, self.parent.name)
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = '12:00:00'
        if 'email' not in kwargs: kwargs['email'] = None
        self.slurm_job = SLURMJob(command, self.p.logs_dir, job_name=str(self), **kwargs)
        self.slurm_job.launch()

    def combine_reads(self):
        shell_output('cat %s > %s' % (' '.join([pool.qiime_fasta.path for pool in self]), self.qiime_reads.path))

    def run_denovo(self): self.otu_denovo.run()
    def run_open_ref(self): self.otu_open_ref.run()
    def run_progressive(self): self.otu_step.run()
    def run_subsample(self): self.otu_subsample.run()
