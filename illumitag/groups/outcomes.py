# Built-in modules #
import os
from collections import Counter, defaultdict, OrderedDict
from commands import getstatusoutput

# Internal modules #
from assemble import Assembled, Unassembled
from common import property_cached, AutoPaths
from fasta.paired import PairedFASTQ

# Third party modules #
import sh

###############################################################################
class BarcodeGroup(PairedFASTQ):
    """A bunch of sequences all having the same type of barcode outcome"""

    all_paths = """
    /fwd.fastq
    /rev.fastq
    /assembled/
    /unassembled/
    """

    def __iter__(self): return iter(self.children)
    def __repr__(self): return '<%s object of pool %i>' % (self.__class__.__name__, self.pool.num)

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        # Paths #
        self.base_dir = self.pool.p.groups_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Super #
        self.fwd_path = self.p.fwd_fastq
        self.rev_path = self.p.rev_fastq
        # Add assembly files #
        self.assembled = Assembled(self)
        self.unassembled = Unassembled(self)
        self.children = (self.assembled, self.unassembled)
        self.first = self.assembled
        # Extra #
        self.samples = self.pool.samples

    def assemble(self):
        command = 'pandaseq -f %s -r %s -u %s -F 1> %s 2> %s'
        command = command % (self.p.fwd_fastq, self.p.rev_fastq, self.unassembled.path, self.assembled.path, self.assembled.p.out)
        getstatusoutput(command)

    def barcode_fastqc(self):
        sh.fastqc(self.fwd_path, '-q')
        os.remove(os.path.splitext(self.fwd_path)[0] + '_fastqc.zip')
        sh.fastqc(self.rev_path, '-q')
        os.remove(os.path.splitext(self.rev_path)[0] + '_fastqc.zip')

    def assembly_fastqc(self):
        sh.fastqc(self.assembled.path, '-q')
        os.remove(os.path.splitext(self.assembled.path)[0] + '_fastqc.zip')

    def check_noalign_counts(self):
        assert self.assembled.stats['noalign'] == self.unassembled.count

###############################################################################
class NoBarcode(BarcodeGroup):
    short_name = "no_barcodes"
    @property
    def doc(self): return "No barcodes found at all"

    @property_cached
    def counter(self):
        return Counter()

#-----------------------------------------------------------------------------#
class OneBarcode(BarcodeGroup):
    short_name = "one_barcodes"
    doc = "Only one barcode found"

    @property_cached
    def counter(self):
        return Counter(p.matches[0].name for p in self.parse_barcodes())

#-----------------------------------------------------------------------------#
class SameBarcode(BarcodeGroup):
    """When the barcodes have a set discrepancy"""
    short_name = "same_barcodes"
    doc = "Both barcodes in same set"

    @property_cached
    def set_counter(self):
        pairs = ((p.matches[0].name, p.matches[1].name) for p in self.parse_barcodes())
        return Counter(map(frozenset, pairs))

    @property_cached
    def counter(self):
        return Counter((m.name for pair in self.parse_barcodes() for m in pair.matches))

#-----------------------------------------------------------------------------#
class BadBarcode(BarcodeGroup):
    """When the barcodes have an index discrepancy"""
    short_name = "bad_barcodes"
    doc = "The two barcodes missmatch"

    @property_cached
    def set_counter(self):
        pairs = ((p.matches[0].name, p.matches[1].name) for p in self.parse_barcodes())
        return Counter(map(frozenset, pairs))

    @property_cached
    def counter(self):
        return Counter((m.name for pair in self.parse_barcodes() for m in pair.matches))

    @property_cached
    def distribution(self):
        result = defaultdict(Counter)
        for k,v in self.set_counter.items():
            first, second = k
            result[first][second] = v
            result[second][first] = v
        return result

#-----------------------------------------------------------------------------#
class GoodBarcode(BarcodeGroup):
    short_name = "good_barcodes"
    doc = "Success: barcodes match"

    @property_cached
    def counter(self):
        return Counter((m.name for pair in self.parse_barcodes() for m in pair.matches))

    @property_cached
    def breakdown(self):
        return OrderedDict([(name, self.counter[name + 'A']) for name in self.barcodes.names])