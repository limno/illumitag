# Built-in modules #

# Internal modules #
from illumitag.helper.chimeras import UchimeRef, UchimeDenovo
from illumitag.fasta.single import FASTA
from illumitag.common.autopaths import AutoPaths
from illumitag.clustering.taxonomy.rdp import SimpleRdpTaxonomy
from illumitag.clustering.taxonomy.crest import SimpleCrestTaxonomy

# Third party modules #

###############################################################################
class Fractions(object):
    """All size fractions from the quality reads"""

    def __repr__(self): return '<%s object of "%s">' % (self.__class__.__name__, self.parent)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)

    all_paths = """
    /low/
    /med/
    /big/
    """

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        # Auto paths #
        self.base_dir = self.parent.p.fractions_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Size fractions #
        self.low = Fraction(self, self.p.low_dir, 430, 446)
        self.med = Fraction(self, self.p.med_dir, 447, 464)
        self.big = Fraction(self, self.p.big_dir, 465, 488)
        self.children = [self.low, self.med, self.big]

###############################################################################
class Fraction(object):
    """One size fraction from the quality reads."""

    all_paths = """
    /reads.fasta
    /refere/
    /denovo/
    /rdp/
    /crest/
    """

    def __init__(self, parent, base_dir, lower_bound, upper_bound):
        # Save parent #
        self.parent, self.fractions = parent, parent
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Bounds #
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        # Size fractions #
        self.reads = FASTA(self.p.reads_fasta)
        self.refere = UchimeRef(self.p.reads, self.p.refere_dir, self)
        self.denovo = UchimeDenovo(self.p.reads, self.p.denovo_dir, self)
        # Classification #
        self.rdp = SimpleRdpTaxonomy(self.reads, self.p.rdp_dir)
        self.crest = SimpleCrestTaxonomy(self.reads, self.p.crest_dir)

    def extract(self):
        self.fractions.pool.quality_reads.only_used.extract_length(self.lower_bound, self.upper_bound, self.reads)
        return self.reads

    def check_chimeras(self):
        self.refere.check()
        self.denovo.check()
        return self.refere, self.denovo

    def check_classification(self):
        self.rdp.assign()
