# Built-in modules #
import sys, os

# Internal modules #
from illumitag.common import moving_average
from illumitag.helper.chimeras import UchimeRef, UchimeDenovo
from illumitag.fasta.single import FASTQ, FASTA
from illumitag.common.autopaths import AutoPaths
from illumitag.common.color import Color

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'
chimera_ref_path = home + 'glob/16s/microbiomeutil-r20110519.fasta'

###############################################################################
class PrimerGroup(object):
    """A bunch of sequences all having the same type of primer outcome
    (and assembly outcome)"""

    all_paths = """
    /orig.fastq
    /n_filtered.fastq
    /qual_filtered.fastq
    /len_filtered.fastq
    /trimmed_barcodes.fasta
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.orig_reads)

    def create(self): self.orig_reads.create()
    def add_read(self, read): self.orig_reads.add_read(read)
    def close(self): self.orig_reads.close()

    def __init__(self, parent):
        # Save parent #
        self.parent, self.assemble_group = parent, parent
        self.samples = parent.samples
        self.pool = self.parent.pool
        # Auto paths #
        self.base_dir = parent.p.groups_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # More #
        self.orig_reads = self.parent.cls(self.p.orig_fastq, samples=self.samples)
        self.n_filtered = self.parent.cls(self.p.n_filtered, samples=self.samples)
        # Quality filtered #
        if self.parent == 'assembled':
            self.qual_filtered = FASTQ(self.p.qual_filtered, samples=self.samples)
            self.len_filtered = FASTQ(self.p.len_filtered_fastq, samples=self.samples)
            self.trimmed_barcodes = FASTA(self.p.trimmed_barcodes)
        # Further #
        self.load()

    def load(self): pass

    def n_filter(self):
        """Called from AssembleGroup.discard_reads_with_n"""
        def no_n_iterator(reads):
            for read in reads:
                if 'N' in read: continue
                yield read
        self.n_filtered.write(no_n_iterator(self.orig_reads))

    def qual_filter(self):
        """Called from Assemble.quality_filter"""
        def good_qual_iterator(reads, threshold=5, windowsize=10):
            for read in reads:
                averaged = moving_average(read.letter_annotations["phred_quality"], windowsize)
                if any([value < threshold for value in averaged]): continue
                yield read
        self.qual_filtered.write(good_qual_iterator(self.n_filtered))

    def len_filter(self):
        """Called from Assemble.length_filter"""
        def good_len_iterator(reads, min_length=300):
            for read in reads:
                if len(read) < min_length: continue
                yield read
        self.len_filtered.write(good_len_iterator(self.qual_filtered))

    def trim_bc(self):
        """Called from Assemble.trim_barcodes"""
        def no_barcodes_iterator(reads):
            for read in reads:
                yield read[self.pool.bar_len:-self.pool.bar_len]
        if self.pool.bar_len == 0:
            self.len_filtered.to_fasta(self.trimmed_barcodes)
        else:
            self.trimmed_barcodes.write(no_barcodes_iterator(self.len_filtered))

###############################################################################
class GoodPrimers(PrimerGroup):
    """Both primers found at right place"""
    short_name = 'good_primers'

    all_paths = PrimerGroup.all_paths + """
    /
    /chimeras_ref/
    /chimeras_denovo/
    /size_fractions/low/low.fasta
    /size_fractions/low/refere/
    /size_fractions/low/denovo/
    /size_fractions/med/med.fasta
    /size_fractions/med/refere/
    /size_fractions/med/denovo/
    /size_fractions/big/big.fasta
    /size_fractions/big/refere/
    /size_fractions/big/denovo/
    """

    def load(self):
        if self.parent == 'assembled':
            # Chimeras #
            self.uchime_ref = UchimeRef(self.trimmed_barcodes.path, self.p.chimeras_ref_dir, self)
            self.uchime_denovo = UchimeDenovo(self.trimmed_barcodes.path, self.p.chimeras_denovo_dir, self)
            # Size fractions #
            self.low_length = FASTA(self.p.low_fasta)
            self.low_refere = UchimeRef(self.p.low_fasta, self.p.low_refere_dir, self)
            self.low_denovo = UchimeDenovo(self.p.low_fasta, self.p.low_denovo_dir, self)
            self.med_length = FASTA(self.p.med_fasta)
            self.med_refere = UchimeRef(self.p.med_fasta, self.p.med_refere_dir, self)
            self.med_denovo = UchimeDenovo(self.p.med_fasta, self.p.med_denovo_dir, self)
            self.big_length = FASTA(self.p.big_fasta)
            self.big_refere = UchimeRef(self.p.big_fasta, self.p.big_refere_dir, self)
            self.big_denovo = UchimeDenovo(self.p.big_fasta, self.p.big_denovo_dir, self)

    def check_chimeras(self):
        # Only assembled sequences #
        if self.parent != 'assembled': return
        # Message #
        message = "----> Checking pool %i %s for chimeras"
        print Color.l_ylw + message % (self.parent.parent.pool.num, self.parent.parent.short_name) + Color.end
        sys.stdout.flush()
        # Check empty #
        if not self.trimmed_barcodes:
            print Color.l_ylw + 'No chimeras to check for %s (empty set)' % (self,) + Color.end
            sys.stdout.flush()
            return
        # Call #
        self.uchime_ref.check()
        self.uchime_denovo.check()

    def check_size_fractions(self):
        self.len_filtered.extract_length(430, 446, self.low_length)
        self.len_filtered.extract_length(447, 464, self.med_length)
        self.len_filtered.extract_length(465, 488, self.big_length)

#-----------------------------------------------------------------------------#
class WrongPrimers(PrimerGroup):
    """Both found but at the wrong place"""
    short_name = 'wrong_primers'

#-----------------------------------------------------------------------------#
class OnlyFwdPrimers(PrimerGroup):
    """Only the forward primer seen"""
    short_name = 'only_fwd_primers'

#-----------------------------------------------------------------------------#
class OnlyRevPrimers(PrimerGroup):
    """Only the reverse primer seen"""
    short_name = 'only_rev_primers'

#-----------------------------------------------------------------------------#
class NoPrimers(PrimerGroup):
    """No primers found"""
    short_name = 'no_primers'
