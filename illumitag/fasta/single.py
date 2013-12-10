# Futures #
from __future__ import division

# Built-in modules #
import os, gzip, re
from collections import Counter, OrderedDict

# Internal modules #
from illumitag.common import isubsample
from illumitag.common.autopaths import FilePath
from illumitag.common.tmpstuff import new_temp_path
from illumitag.helper.barcodes import ReadWithBarcodes
from illumitag.helper.primers import ReadWithPrimers
from illumitag.common.cache import property_cached
from illumitag.common.color import Color

# Third party modules #
import sh, shutil
from Bio import SeqIO

################################################################################
class FASTA(FilePath):
    """A single FASTA file somewhere in the filesystem"""
    extension = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0

    def __init__(self, path, samples=None, primers=None):
        # Basic #
        self.path = path
        # Extra #
        self.samples = samples
        self.primers = primers

    @property
    def first_read(self):
        self.open()
        seq = SeqIO.parse(self.handle, self.extension).next()
        self.close()
        return seq

    @property_cached
    def count(self):
        """Should probably check file size instead of just caching once #TODO"""
        if self.path.endswith('gz'): return int(sh.zgrep('-c', "^>", self.path, _ok_code=[0,1]))
        else: return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    def open(self):
        if self.path.endswith('gz'): self.handle = gzip.open(self.path, 'r')
        else: self.handle = open(self.path, 'r')

    def close(self):
        if hasattr(self, 'buffer'): self.flush()
        self.handle.close()

    def copy(self, path):
        shutil.copy2(self.path, path)

    def create(self):
        self.buffer = []
        self.buf_count = 0
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        self.handle = open(self.path, 'w')

    def add_read(self, read):
        self.buffer.append(read)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def add_seqrecord(self, seqrecord):
        self.buffer.append(seqrecord)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def flush(self):
        for read in self.buffer: SeqIO.write(read, self.handle, self.extension)
        self.buffer = []

    def write(self, reads):
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        with open(self.path, 'w') as self.handle: SeqIO.write(reads, self.handle, self.extension)

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.extension)

    def parse_barcodes(self):
        return (ReadWithBarcodes(r, self.samples) for r in self.parse())

    def parse_primers(self):
        return (ReadWithPrimers(r, self.primers) for r in self.parse())

    @property_cached
    def barcode_counter(self):
        return Counter((str(m) for read in self.parse_barcodes() for m in read.matches))

    @property_cached
    def good_barcodes_breakdown(self):
        return OrderedDict([(name, self.barcode_counter[name + 'F']) for name in self.samples.bar_names])

    @property_cached
    def lengths(self):
        return Counter((len(s) for s in self.parse()))

    def shorter_than(self, value):
        return 100 * sum((v for k,v in self.lengths.items() if k < value)) / self.count

    def subsample(self, down_to, new_path=None):
        # Auto path #
        if not new_path: new_path = self.p.subsample
        # Check size #
        if down_to > len(self):
            message = "Can't subsample %s down to %i. Only down to %i."
            print Color.ylw + message % (self, down_to, len(self)) + Color.end
            self.copy(new_path)
            return
        # Make new file #
        self.subsampled = self.__class__(new_path)
        self.subsampled.create()
        # Do it #
        for seq in isubsample(self, down_to):
            self.subsampled.add_seqrecord(seq)
        # Clean up #
        self.subsampled.close()
        # Did it work #
        assert len(self.subsampled) == down_to

    def rename_with_num(self, prefix, new_path=None, remove_desc=True):
        # Temporary path #
        if new_path is None: numbered = self.__class__(new_temp_path())
        else: numbered = self.__class__(new_path)
        # Generator #
        def numbered_iterator():
            for i,read in enumerate(self):
                read.id = prefix + str(i)
                if remove_desc: read.description = ""
                yield read
        # Do it #
        numbered.write(numbered_iterator())
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            shutil.move(numbered, self.path)

#-----------------------------------------------------------------------------#
class FASTQ(FASTA):
    """A single FASTQ file somewhere in the filesystem"""
    extension = 'fastq'

    @property_cached
    def count(self):
        if self.path.endswith('gz'): return int(sh.zgrep('-c', "^+$", self.path, _ok_code=[0,1]))
        return int(sh.grep('-c', "^+$", self.path, _ok_code=[0,1]))

    def to_fasta(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'fasta')

    def fastqc(self):
        # Call #
        sh.fastqc(self.path, '-q')
        # Paths #
        zip_file = self.prefix_path + '_fastqc.zip'
        report_dir = self.prefix_path + '_fastqc/'
        #images_dir = self.prefix_path + '_fastqc/Images/'
        # Clean up #
        os.remove(zip_file)
        # Return #
        return report_dir

#-----------------------------------------------------------------------------#
class SizesFASTA(FASTA):
    """A Fasta with cluster weights"""

    @property_cached
    def count(self):
        get_size = lambda x: int(re.findall("size=([0-9]+)", x)[0])
        sizes = (get_size(r.description) for r in self)
        return sum(sizes)