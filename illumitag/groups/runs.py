# Built-in modules #

# Third party modules #

# Internal modules #

###############################################################################
class Run(object):
    """An illumina run containing several pools."""

    all_paths = """
    /
    """

    def __repr__(self): return '<%s object number %i>' % (self.__class__.__name__, self.num)
    def __iter__(self): return iter(self.pools)
    def __len__(self): return self.count
    def __getitem__(self, key): return self.samples[key]

    @property
    def first(self): return self.pools[0]

    @property
    def label(self): return self.first.run_label

    def __init__(self, num, pools):
        self.num = num
        self.pools = pools

    def load(self):
        for p in self.pools: p.load()