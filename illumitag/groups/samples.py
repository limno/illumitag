# Built-in modules #

# Internal modules #
from illumitag.fasta.single import FASTQ, FASTA
from illumitag.helper.barcodes import bar_len
from illumitag.common.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Samples(object):
    """A list-like container for all samples of a pool"""

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)
    def __getitem__(self, key): return self.children[key]

    @property
    def first(self): return self.children[0]

    def __init__(self, parent):
        # Basic #
        self.parent, self.pool = parent, parent
        self.info = parent.info['samples']
        self.children = [Sample(s, self) for s in self.info]
        # Sort them #
        self.children.sort(key=lambda x: x.num)
        # Check number are all there #
        assert [s.num for s in self] == range(1,51)
        # Check names are unique #
        names = [s.short_name for s in self if s.used]
        assert len(names) == len(set(names))

    def load(self):
        # Barcodes #
        self.bars_F = [s.fwd_str for s in self]
        self.bars_R = [s.rev_str for s in self]
        # Barcode names #
        self.bar_names = ['barcode%i' % s.num for s in self]
        self.bar_sided_names = [name + side for name in self.bar_names for side in ('F','R')]
        self.bar_names_F = [name + 'F' for name in self.bar_names]
        self.bar_names_R = [name + 'R' for name in self.bar_names]
        self.all_bar_pairs = [(a,b) for a in self.bar_sided_names for b in self.bar_sided_names if a[:-1] != b[:-1]]
        # Primer size #
        self.trim_fwd = bar_len + self.pool.primers.fwd_len
        self.trim_rev = bar_len + self.pool.primers.rev_len
        # Children #
        for s in self: s.load()

###############################################################################
class Sample(FASTQ):
    """All sequences with the same barcode pair grouped together"""

    all_paths = """
    /orig.fastq
    /trimmed.fastq
    /renamed.fastq
    /reads.fasta
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __str__(self): return self.bar_name

    def __init__(self, info, parent):
        # Save attributes #
        self.info = info
        self.parent = parent
        self.pool = parent.pool
        # Basic #
        self.short_name = info['name']
        self.group_name = info['group']
        self.num = int(info['num'])
        self.used = bool(info['used'])
        self.fwd_str = info['fwd']
        self.rev_str = info['rev']
        # Other #
        self.bar_name = 'barcode%i' % self.num
        self.name = 'run%i_pool%i_sample%i' % (self.pool.run_num, self.pool.num, self.num)

    def load(self):
        # Paths #
        self.base_dir = self.pool.p.samples_dir + self.bar_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.path = str(self.p.orig_fastq)
        # Distances #
        self.trim_fwd = self.pool.samples.trim_fwd
        self.trim_rev = self.pool.samples.trim_rev
        # Files #
        self.trimmed = FASTQ(self.p.trimmed)
        self.renamed = FASTQ(self.p.renamed)
        self.fasta = FASTA(self.p.reads_fasta)

    def process(self):
        def no_primers_iterator(reads):
            for read in reads:
                yield read[self.trim_fwd:-self.trim_rev]
        self.trimmed.write(no_primers_iterator(self))
        self.trimmed.rename_with_num(self.name + '_read', self.renamed)
        self.renamed.to_fasta(self.fasta)
