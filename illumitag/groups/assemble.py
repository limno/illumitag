# Futures #
from __future__ import division

# Built-in modules #
import re
from collections import Counter

# Internal modules #
from primers import GoodPrimers, WrongPrimers, OnlyFwdPrimers, OnlyRevPrimers, NoPrimers
from illumitag.common import tail, flatten, reverse_compl_with_name
from illumitag.fasta.single import FASTQ, FASTA
from illumitag.common.cache import property_cached
from illumitag.common.autopaths import AutoPaths

# Third party modules #

###############################################################################
class AssembleGroup(object):
    """A bunch of sequences all having the same type of assembly outcome
    (and barcode outcome)"""

    def __iter__(self): return iter(self.children)
    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return self.count
    def __ne__(self, other): return not self.__eq__(other)

    def __init__(self, parent):
        # Save parent #
        self.parent, self.outcome = parent, parent
        self.samples = parent.samples
        # Extra #
        self.pool = self.outcome.parent
        self.samples = self.pool.samples
        self.primers = self.pool.primers
        # Load #
        self.load()
        # Auto paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # All primer outcomes #
        self.good_primers     = GoodPrimers(self)
        self.wrong_primers    = WrongPrimers(self)
        self.only_fwd_primers = OnlyFwdPrimers(self)
        self.only_rev_primers = OnlyRevPrimers(self)
        self.no_primers       = NoPrimers(self)
        # Group them #
        self.children = (self.good_primers, self.wrong_primers, self.only_fwd_primers, self.only_rev_primers, self.no_primers)
        self.first = self.good_primers

    @property
    def flipped_iterator(self):
        for r in self.parse_barcodes():
            if r.first.set == 'R' or r.last.set == 'F': yield reverse_compl_with_name(r.read)
            else: yield r.read

    def flip_reads(self):
        self.flipped_reads.write(self.flipped_iterator)

    def dont_flip_reads(self):
        self.flipped_reads.link_from(self.path)

    def make_primer_groups(self):
        bar_len = self.parent.parent.bar_len
        for g in self.children: g.create()
        for r in self.flipped_reads.parse_primers():
            if r.fwd_pos is not None and r.rev_pos is not None:
                if r.fwd_pos == bar_len and r.rev_pos == -bar_len: self.good_primers.add_read(r.read)
                else:                                              self.wrong_primers.add_read(r.read)
            elif r.fwd_pos:                                        self.only_fwd_primers.add_read(r.read)
            elif r.rev_pos:                                        self.only_rev_primers.add_read(r.read)
            else:                                                  self.no_primers.add_read(r.read)
        for g in self.children: g.close()

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

    def load(self):
        self.cls = FASTQ
        self.base_dir = self.outcome.p.assembled_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.path = self.p.orig_fastq
        self.flipped_reads = FASTQ(self.p.flipped, self.samples, self.primers)

    @property_cached
    def stats(self):
        result = {}
        result['raw'] = tail(self.p.out)
        if "pandaseq: error" in result['raw']: raise Exception("Pandaseq did not run properly")
        result['distrib'] = re.findall('STAT\tOVERLAPS\t(.+)$', result['raw'], re.M)
        result['distrib'] = map(int, result['distrib'][0].split())
        result['lengths'] = flatten([[i+1]*v for i,v in enumerate(result['distrib'])])
        result['noalign'] = int(re.findall('STAT\tNOALGN\t(.+)$', result['raw'], re.M)[0])
        result['lowqual'] = int(re.findall('STAT\tLOWQ\t(.+)$', result['raw'], re.M)[0])
        result['loss'] = 100 * sum(result['distrib'][100:]) / sum(result['distrib'])
        return result

    def quality_filter(self):
        for g in self.children: g.qual_filter()

    def length_filter(self):
        for g in self.children: g.len_filter()

    def trim_barcodes(self):
        for g in self.children: g.trim_bc()

#-----------------------------------------------------------------------------#
class Unassembled(AssembleGroup, FASTA):
    def __eq__(self, other): return other == 'unassembled'

    all_paths = """
    /orig.fasta
    /flipped.fasta
    /groups/
    """

    def load(self):
        self.cls = FASTA
        self.base_dir = self.outcome.p.unassembled_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.path = self.p.orig_fasta
        self.flipped_reads = FASTA(self.p.flipped, self.samples, self.primers)
