# Built-in modules #
import shutil

# Internal modules #
import illumitag
from illumitag.common.tmpstuff import new_temp_path
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.running.cluster_runner import ClusterRunner
from illumitag.clustering.otu.uparse import UparseOTUs

# Third party modules #
import pandas
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

    def __init__(self, samples, name, base_dir=None):
        # Save samples #
        self.name = name
        self.samples, self.children = samples, samples
        # Check names are unique #
        names = [s.short_name for s in samples if s.used]
        assert len(names) == len(set(names))
        # Figure out pools #
        self.pools = list(set([s.pool for s in self.samples]))
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

    def run(self, *args, **kwargs):
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        self.runner.run_slurm()

    def combine_reads(self):
        paths = [sample.fasta.path for sample in self]
        shell_output('cat %s > %s' % (' '.join(paths), self.reads))

    def set_size(self, length=None):
        if length is None: return
        self.size_set = FASTA(new_temp_path())
        def sized_iterator(reads):
            for read in reads:
                if len(read) < length: continue
                yield read[:length]
        self.size_set.write(sized_iterator(self.reads))
        self.size_set.close()
        # Replace it #
        self.reads.remove()
        shutil.move(self.size_set, self.reads)

    def run_uparse(self): self.otu_uparse.run()

    @property
    def metadata(self):
        return pandas.DataFrame([s.info for s in self], index=[s.short_name for s in self])

    def export_metadata(self):
        self.metadata.to_csv(self.p.metadata, sep='\t', encoding='utf-8')


