# Built-in modules #
import re
from collections import Counter

# Internal modules #
from primers import GoodPrimers, WrongPrimers, OnlyFwdPrimers, OnlyRevPrimers, NoPrimers
from illumitag.common import property_cached, AutoPaths, tail, flatten, reverse_compl_with_name
from illumitag.fasta.single import FASTQ, FASTA
from illumitag.helper.barcodes import bar_len

# Third party modules #
from Bio import SeqIO

###############################################################################
class AssembleGroup(object):
    """A bunch of sequences all having the same type of assembly outcome
    (and barcode outcome)"""

    def __iter__(self): return iter(self.children)
    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return self.count

    def load(self):
        # All primer outcomes #
        self.good_primers     = GoodPrimers(self)
        self.wrong_primers    = WrongPrimers(self)
        self.only_fwd_primers = OnlyFwdPrimers(self)
        self.only_rev_primers = OnlyRevPrimers(self)
        self.no_primers       = NoPrimers(self)
        self.children = (self.good_primers, self.wrong_primers, self.only_fwd_primers, self.only_rev_primers, self.no_primers)
        self.first = self.good_primers
        # Extra #
        self.pool = self.outcome.pool
        self.samples = self.pool.samples
        self.primers = self.pool.primers

    def make_primer_groups(self):
        for g in self.children: g.create()
        for r in self.flipped_reads.parse_primers():
            if r.fwd_pos and r.rev_pos:
                if r.fwd_pos == bar_len and r.rev_pos == -bar_len: self.good_primers.add_read(r.read)
                else:                                              self.wrong_primers.add_read(r.read)
            elif r.fwd_pos:                                        self.only_fwd_primers.add_read(r.read)
            elif r.rev_pos:                                        self.only_rev_primers.add_read(r.read)
            else:                                                  self.no_primers.add_read(r.read)
        for g in self.children: g.close()

    @property
    def flipped_iterator(self):
        for r in self.parse_barcodes():
            if r.first.set == 'R' or r.last.set == 'F': yield reverse_compl_with_name(r.read)
            else: yield r.read

    @property_cached
    def primer_positions(self):
        # Count positions #
        all_fwd_pos, all_rev_pos = Counter(), Counter()
        for r in self.flipped_reads.parse_primers():
            if r.fwd_pos: all_fwd_pos.update((r.fwd_pos,))
            if r.rev_pos: all_rev_pos.update((r.rev_pos,))
        # Return results #
        return all_fwd_pos, all_rev_pos

    def discard_reads_with_n(self):
        for g in self.children: g.n_filter()

###############################################################################
class Assembled(AssembleGroup, FASTQ):
    all_paths = """
    /orig.fastq
    /flipped.fastq
    /pandaseq.out
    /groups/
    """

    def __eq__(self, other): return other == 'assembled'
    def __ne__(self, other): return not self.__eq__(other)

    def __init__(self, parent):
        # Save parent #
        self.parent, self.outcome = parent, parent
        self.samples = parent.samples
        # Auto paths #
        self.base_dir = self.outcome.p.assembled_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Super #
        self.load()
        # Make fastq files #
        self.path = self.p.orig_fastq
        self.flipped_reads = FASTQ(self.p.flipped, self.samples, self.primers)

    @property_cached
    def stats(self):
        result = {}
        result['raw'] = tail(self.p.out)
        result['distrib'] = re.findall('^thread[0-99]\tSTAT\tOVERLAPS\t(.+)$', result['raw'], re.M)
        result['distrib'] = map(int, result['distrib'][0].split())
        result['lengths'] = flatten([[i+1]*v for i,v in enumerate(result['distrib'])])
        result['noalign'] = int(re.findall('\tSTAT\tNOALGN\t(.+)$', result['raw'], re.M)[0])
        return result

    def flip_reads(self):
        with open(self.flipped_reads.path, 'w') as handle: SeqIO.write(self.flipped_iterator, handle, 'fastq')

    def quality_filter(self):
        for g in self.children: g.qual_filter()

    def len_filter(self):
        for g in self.children: g.len_filter()

#-----------------------------------------------------------------------------#
class Unassembled(AssembleGroup, FASTA):
    all_paths = """
    /orig.fasta
    /flipped.fasta
    /groups/
    """

    def __eq__(self, other): return other == 'unassembled'
    def __ne__(self, other): return not self.__eq__(other)

    def __init__(self, parent):
        # Save parent #
        self.parent, self.outcome = parent, parent
        self.samples = parent.samples
        # Auto paths #
        self.base_dir = self.outcome.p.unassembled_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Super #
        self.load()
        # Make fasta files #
        self.path = self.p.orig_fasta
        self.flipped_reads = FASTA(self.p.flipped, self.samples, self.primers)

    def flip_reads(self):
        with open(self.flipped_reads.path, 'w') as handle: SeqIO.write(self.flipped_iterator, handle, 'fasta')
