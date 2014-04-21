# Built-in modules #
import shutil

# Internal modules #
import illumitag
from illumitag.common.tmpstuff import new_temp_path
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.running.cluster_runner import ClusterRunner
from illumitag.clustering.otu.uparse import UparseOTUs
from illumitag.clustering.otu.uclust import UclustOTUs
from illumitag.clustering.reporting import ClusterReporter

# Third party modules #
import pandas
from tqdm import tqdm
from shell_command import shell_output

###############################################################################
class Cluster(object):
    """Analyzes a group of samples."""

    all_paths = """
    /reads/all_reads.fasta
    /otus/
    /logs/
    /metadata.csv
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % (self.__class__.__name__, self.name, len(self.samples))
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key): return self.samples[key]

    @property
    def first(self): return self.children[0]

    @property
    def count_seq(self):
        return sum([len(sample) for sample in self])

    def __init__(self, samples, name, base_dir=None):
        # Save samples #
        self.name = name
        self.samples, self.children = samples, samples
        # Check names are unique #
        names = [s.short_name for s in samples if s.used]
        assert len(names) == len(set(names))
        # Figure out pools #
        self.pools = list(set([s.pool for s in self.samples]))
        self.pools.sort(key = lambda x: x.id_name)
        # Load them #
        for p in self.pools: p.load()
        # Dir #
        if base_dir: self.base_dir = base_dir
        else: self.base_dir = illumitag.view_dir + "clusters/" + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner = ClusterRunner(self)
        # FASTA #
        self.reads = FASTA(self.p.all_reads_fasta)
        # OTU picking #
        self.otu_uparse = UparseOTUs(self)
        self.otu_uclust = UclustOTUs(self)
        # Reporting #
        self.reporter = ClusterReporter(self)

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        self.runner.run_slurm()

    def process_samples(self):
        for sample in tqdm(self): sample.process()

    def combine_reads(self):
        paths = [sample.fasta.path for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))

    def set_size(self, length):
        """Trim all sequences to a specific length starting from the end"""
        self.size_trimmed = FASTA(new_temp_path())
        def trim_iterator(reads):
            for read in reads:
                if len(read) < length: continue
                yield read[-length:]
        self.size_trimmed.write(trim_iterator(self.reads))
        self.size_trimmed.close()
        # Replace it #
        self.reads.remove()
        shutil.move(self.size_trimmed, self.reads)

    def run_uparse(self): self.otu_uparse.run()

    @property
    def metadata(self):
        return pandas.DataFrame([s.info for s in self], index=[s.short_name for s in self])

    def export_metadata(self):
        self.metadata.to_csv(self.p.metadata, sep='\t', encoding='utf-8')