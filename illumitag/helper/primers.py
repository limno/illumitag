# Built-in modules #
import re

# Internal modules #

# Third party modules #
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Constants #
iupac = {'A':'A', 'G':'G', 'T':'T', 'C':'C', 'M':'AC', 'R':'AG', 'W':'AT', 'S':'CG', 'Y':'CT', 'K':'GT', 'V':'ACG', 'H':'ACT', 'D':'AGT', 'B':'CGT', 'X':'ACGT', 'N':'ACGT'}

###############################################################################
class TwoPrimers(object):
    """A container for the two primers of a pool"""

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent.id_name)
    def __len__(self): return 2

    def __init__(self, parent):
        self.parent, self.pool = parent, parent
        self.info = parent.info['primers']
        # Basic #
        self.name = self.info['name']
        self.fwd_name = self.info['forward']['name']
        self.rev_name = self.info['reverse']['name']
        self.fwd_str = ['forward']['sequence']
        self.rev_str = ['reverse']['sequence']

    def load(self):
        self.fwd_seq = Seq(self.fwd_str, IUPAC.ambiguous_dna)
        self.rev_seq = Seq(self.rev_str, IUPAC.ambiguous_dna)
        self.fwd_regex = re.compile(''.join(['[' + iupac[char] + ']' for char in self.fwd_seq]))
        self.rev_regex = re.compile(''.join(['[' + iupac[char] + ']' for char in self.rev_seq.reverse_complement()]))

###############################################################################
class ReadWithPrimers(object):
    def __init__(self, read, primers):
        self.read = read
        self.fwd_match = primers.fwd_regex.search(str(read.seq))
        self.rev_match = primers.rev_regex.search(str(read.seq))
        self.fwd_pos = self.fwd_match.start() if self.fwd_match else None
        self.rev_pos = self.rev_match.end() - len(read) if self.rev_match else None
