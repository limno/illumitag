# Built-in modules #
import os, sys, gzip
from itertools import izip

# Internal modules #
from common import property_cached
from helper.barcodes import ReadPairWithBarcode

# Third party modules #
import sh
from Bio import SeqIO

###############################################################################
class PairedFASTQ(object):
    """Read and write FASTQ file pairs without using too much RAM"""
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, fwd_path, rev_path, samples=None, primers=None):
        # Basic #
        self.fwd_path = fwd_path
        self.rev_path = rev_path
        # Extra #
        self.samples = samples
        self.primers = primers

    @property_cached
    def count(self):
        return int(sh.grep('-c', "^+$", self.fwd_path, _ok_code=[0,1]))

    def open(self):
        # Support for compression #
        if self.fwd_path.endswith('gz'):
            self.fwd_handle = gzip.open(self.fwd_path, 'r')
        else: self.fwd_handle = open(self.fwd_path, 'r')
        if self.rev_path.endswith('gz'):
            self.rev_handle = gzip.open(self.rev_path, 'r')
        else: self.rev_handle = open(self.rev_path, 'r')

    def close(self):
        self.flush()
        self.fwd_handle.close()
        self.rev_handle.close()

    def create(self):
        # The buffer #
        self.buffer = []
        self.buf_count = 0
        # Directory #
        self.fwd_dir = os.path.dirname(self.fwd_path)
        self.rev_dir = os.path.dirname(self.rev_path)
        if not os.path.exists(self.fwd_dir):
            os.makedirs(self.fwd_dir)
        if not os.path.exists(self.rev_dir):
            os.makedirs(self.rev_dir)
        # The files #
        self.fwd_handle = open(self.fwd_path, 'w')
        self.rev_handle = open(self.rev_path, 'w')

    def add_pair(self, pair):
        self.buffer.append(pair)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0:
            sys.stderr.write('.')
            self.flush()

    def flush(self):
        for pair in self.buffer:
            SeqIO.write(pair.fwd, self.fwd_handle, 'fastq')
            SeqIO.write(pair.rev, self.rev_handle, 'fastq')
        self.buffer = []

    def parse(self):
        self.open()
        return izip(SeqIO.parse(self.fwd_handle, 'fastq'),
                    SeqIO.parse(self.rev_handle, 'fastq'))

    def parse_barcodes(self):
        return (ReadPairWithBarcode(f, r, self.samples) for f,r in self.parse())