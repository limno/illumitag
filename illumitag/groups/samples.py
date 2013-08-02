# Built-in modules #

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.fasta.single import FASTQ, FASTA
from illumitag.fasta.other import GroupFile
from illumitag.helper.barcodes import bar_len

# Third party modules #

###############################################################################
class Samples(object):
    """A list-like container for all samples of a pool"""

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent.id_name)
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
        # Children #
        for s in self: s.load()
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

###############################################################################
class Sample(FASTQ):
    """All sequences with the same barcode pair grouped together"""

    all_paths = """
    /reads.fastq
    /subsample.fastq
    /mothur/reads.fasta
    /mothur/groups.tsv
    """

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent.id_name)
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

    def load(self):
        # Paths #
        self.base_dir = self.pool.p.samples_dir + self.bar_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Super #
        self.path = self.p.reads_fastq
        # Subsample #
        self.subsampled = FASTQ(self.p.subsample_fastq)
        # Mothur #
        self.mothur_fasta = FASTA(self.p.mothur_fasta)
        self.mothur_groups = GroupFile(self.p.mothur_groups)

    def make_mothur_output(self):
        # Trimmed fasta #
        self.mothur_fasta.create()
        for r in self.subsampled: self.mothur_fasta.add_seqrecord(r[self.pool.trim_start:-self.pool.trim_end])
        self.mothur_fasta.close()
        # The groups file #
        self.mothur_groups.create()
        group_name = "\tpool%i-%s\n" % (self.pool.num, self.name)
        for r in self.subsampled: self.mothur_groups.handle.write(r.id + group_name)
        self.mothur_groups.close()
