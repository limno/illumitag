# Built-in modules #
import os
from collections import defaultdict

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.fasta.other import QualFile, GroupFile
from illumitag.helper.barcodes import bar_len

# Third party modules #
from Bio.SeqIO.FastaIO import FastaWriter

# Constants #

###############################################################################
class QualityReads(object):
    """A set of sequences determined to be quality controlled"""

    all_paths = """
    /mothur_reads.fasta
    /mothur_reads.qual
    /mothur_groups.tsv
    /qiime_reads.fasta
    /only_used_samples.fasta
    /trimmed.fasta
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.orig_reads)

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        self.samples = parent.samples
        # Auto paths #
        self.base_dir = parent.p.quality_dir + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        self.untrimmed = self.pool.good_barcodes.assembled.good_primers.len_filtered
        self.only_used = FASTA(self.p.only_used, samples=self.samples)
        self.trimmed = FASTA(self.p.trimmed)
        # Qiime output #
        self.qiime_fasta = FASTA(self.p.qiime_fasta)
        # Mothur #
        self.mothur_fasta = FASTA(self.p.mothur_fasta)
        self.mothur_qual = QualFile(self.p.mothur_qual)
        self.mothur_groups = GroupFile(self.p.mothur_groups)
        # Primer size #
        self.trim_fwd = self.pool.samples.trim_fwd
        self.trim_rev = self.pool.samples.trim_rev

    def filter_unused(self):
        def no_unused_iterator(reads):
            for r in reads.parse_barcodes():
                if r.first.sample.used: yield r.read
        self.only_used.write(no_unused_iterator(self.untrimmed))

    def trim_primers(self):
        def no_primers_iterator(reads):
            for read in reads:
                yield read[self.trim_fwd:-self.trim_rev]
        self.trimmed.write(no_primers_iterator(self.only_used))

    def make_mothur_output(self):
        # Trimmed fasta #
        if os.path.exists(self.mothur_fasta.path): os.remove(self.mothur_fasta.path)
        os.symlink(self.trimmed.path, self.mothur_fasta.path)
        # The groups file #
        self.mothur_groups.create()
        for r in self.untrimmed.parse_barcodes():
            sample_name = r.first.sample.short_name
            read_name = '%s\t%s\n' % (r.read.id, sample_name)
            self.mothur_groups.handle.write(read_name)
        self.mothur_groups.close()

    def make_qiime_output(self):
        # Prepare fasta writer #
        handle = open(self.qiime_fasta.path, 'w')
        writer = FastaWriter(handle, wrap=0)
        writer.write_header()
        # Counter #
        counter = defaultdict(int)
        # Do it #
        for r in self.only_used.parse_barcodes():
            sample_name = r.first.sample.short_name
            counter[sample_name] += 1
            r.read.id = '%s_%i %s' % (sample_name, counter[sample_name], r.read.id)
            bar_seq = r.read.seq[0:bar_len]
            r.read.description = "orig_bc=%s new_bc=%s bc_diffs=0" % (bar_seq, bar_seq)
            writer.write_record(r.read[self.trim_fwd:-self.trim_rev])
        # Close #
        writer.write_footer()
        handle.close()