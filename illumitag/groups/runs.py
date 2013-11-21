# Built-in modules #

# Internal modules #
from aggregate import Collection, Aggregate
from illumitag.common.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Runs(Collection):
    """A collection of runs."""
    pass

###############################################################################
class Run(Aggregate):
    """An illumina run containing several pools."""

    def __repr__(self): return '<%s object number %i with %i pools>' % \
                               (self.__class__.__name__, self.num, len(self))

    @property
    def label(self): return self.first.run_label

    def __init__(self, num, pools, out_dir):
        # Attributes #
        self.num = num
        self.name = "run%i" % num
        self.pools, self.children = pools, pools
        self.loaded = False
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra #
        self.meta_data_path = None