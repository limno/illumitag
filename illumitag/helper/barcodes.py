# Built-in modules #

# Third party modules #

# Constants #
bar_len = 7

###############################################################################
class BarcodeMatch(object):
    """Given a 7 nucleotide sequence and a collection of samples,
    will find the match if it exists"""

    def __nonzero__(self): return bool(self.set)
    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, str(self))
    def __str__(self): return str(self.sample) + self.set

    def __init__(self, bar, samples):
        # Attributes #
        self.bar = bar
        # Default values #
        index_F, index_R = -1, -1
        self.set, self.sample = None, None
        # Search #
        try: index_F = samples.bars_F.index(bar)
        except ValueError: pass
        try: index_R = samples.bars_R.index(bar)
        except ValueError: pass
        # Record #
        if index_F is not -1:
            self.set = "F"
            self.sample = samples[index_F]
        if index_R is not -1:
            self.set = "R"
            self.sample = samples[index_R]

###############################################################################
class ReadWithBarcodes(object):
    def __init__(self, read, samples):
        self.read = read
        self.first = BarcodeMatch(read.seq.tostring()[0:bar_len], samples)
        self.last = BarcodeMatch(read.reverse_complement().seq.tostring()[0:bar_len], samples)
        self.matches = (self.first, self.last)

###############################################################################
class ReadPairWithBarcode(object):
    def __init__(self, fwd, rev, samples):
        self.fwd = fwd
        self.rev = rev
        self.samples = samples

    @property
    def matches(self):
        fwd_m = BarcodeMatch(self.fwd.seq.tostring()[0:bar_len], self.samples)
        rev_m = BarcodeMatch(self.rev.seq.tostring()[0:bar_len], self.samples)
        return [m for m in (fwd_m,rev_m) if m]

    @property
    def illumina_mid(self):
        return self.fwd.description[-16:]