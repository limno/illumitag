# Built-in modules #
import sys, os

# Internal modules #
from illumitag.common import AutoPaths, moving_average, Color
from illumitag.helper.chimeras import UchimeRef, UchimeDenovo
from illumitag.fasta.single import FASTQ, FASTA, BarcodedFASTQ

# Third party modules #
import sh
#from shell_command import shell_output

# Constants #
home = os.environ['HOME'] + '/'
chimera_ref_path = home + 'ILLUMITAQ/dataset1/raw/database/silvagold.fasta'
chimera_ref_path = home + 'ILLUMITAQ/combined-qiime/info/microbiomeutil-r20110519.fasta'

###############################################################################
class PrimerGroup(object):
    """A bunch of sequences all having the same type of primer outcome
    (and assembly outcome)"""
    short_name = 'primers_group'

    all_paths = """
    /orig.fastq
    /n_filtered.fastq
    /qual_filtered.fastq
    /len_filtered.fastq
    /trimmed.fastq
    /trimmed.fasta
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.orig_reads)

    def __init__(self, parent):
        # Save parent #
        self.parent = parent
        self.samples = parent.samples
        self.pool = self.parent.pool
        # Auto paths #
        self.base_dir = parent.p.groups_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # More #
        if self.parent == 'assembled': FASTCLASS = FASTQ
        else:                          FASTCLASS = FASTA
        self.orig_reads = FASTCLASS(self.p.orig_fastq, samples=self.samples)
        self.n_filtered = FASTCLASS(self.p.n_filtered, samples=self.samples)
        # Quality filtered #
        if self.parent == 'assembled':
            self.qual_filtered = FASTQ(self.p.qual_filtered, samples=self.samples)
            self.len_filtered = BarcodedFASTQ(self.p.len_filtered_fastq, samples=self.samples)
            self.trimmed_fastq = FASTQ(self.p.trimmed_fastq)
            self.trimmed_fasta = FASTA(self.p.trimmed_fasta)
        # Extra #
        self.load()

    def load(self): pass

    def fastqc(self):
        sh.fastqc(self.orig_reads.path, '-q')
        os.remove(os.path.splitext(self.orig_reads.path)[0] + '_fastqc.zip')

    def n_filter(self):
        def no_n_iterator(reads):
            for read in reads:
                if 'N' in read: continue
                yield read
        self.n_filtered.write(no_n_iterator(self.orig_reads))

    def qual_filter(self):
        def good_qual_iterator(reads, threshold=5, windowsize=10):
            for read in reads:
                averaged = moving_average(read.letter_annotations["phred_quality"], windowsize)
                if any([value < threshold for value in averaged]): continue
                yield read
        self.qual_filtered.write(good_qual_iterator(self.n_filtered))

    def len_filter(self):
        def good_len_iterator(reads, max_overlap=100):
            min_length = 250 + 250 - max_overlap
            for read in reads:
                if len(read) < min_length: continue
                yield read
        self.len_filtered.write(good_len_iterator(self.qual_filtered))

    def trim_barcodes(self):
        def no_bars_iterator(reads):
            for read in reads:
                yield read[self.pool.trim_fwd:-self.pool.trim_rev]
        self.trimmed_fastq.write(no_bars_iterator(self.len_filtered))
        self.trimmed_fastq.to_fasta(self.trimmed_fasta.path)

    def create(self): self.orig_reads.create()
    def add_read(self, read): self.orig_reads.add_read(read)
    def close(self): self.orig_reads.close()

###############################################################################
class GoodPrimers(PrimerGroup):
    """Both primers found at right place"""
    short_name = 'good_primers'

    all_paths = PrimerGroup.all_paths + """
    /
    /chimeras_ref/
    /chimeras_denovo/
    """

    def load(self):
        self.uchime_ref = UchimeRef(self.p.trimmed_fasta, self.p.chimeras_ref_dir, self)
        self.uchime_denovo = UchimeDenovo(self.p.trimmed_fasta, self.p.chimeras_denovo_dir, self)

    def check_chimeras(self):
        # Only assembled sequences #
        if self.parent != 'assembled': return
        # Message #
        message = "----> Checking pool %i %s for chimeras"
        print Color.l_ylw + message % (self.parent.parent.pool.num, self.parent.parent.short_name) + Color.end
        sys.stdout.flush()
        # Check empty #
        if not self.uchime_ref.fasta:
            print Color.l_ylw + 'No chimeras to check for %s (empty set)' % (self,) + Color.end
            sys.stdout.flush()
            return
        # Call #
        self.uchime_ref.check()
        self.uchime_denovo.check()

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
