# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.graphs import aggregate_plots
from illumitag.analysis import Analysis

# Third party modules #
import sh

###############################################################################
class Collection(object):
    """A collection of aggregates."""

    def __repr__(self): return 'Collection: %s' % (self.children)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)

    def __init__(self, children):
        self.children = children

    @property
    def first(self): return self.children[0]

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return [c for c in self.children if c.name == key][0]
        elif isinstance(key, int):
            if hasattr(self.first, 'num'): return [c for c in self.children if c.num == key][0]
            else: return self.children[key]
        else: raise TypeError('key')

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several pools."""

    all_paths = """
    /graphs/
    """

    def __repr__(self): return '<%s object "%s" with %i pools>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.pools)
    def __len__(self): return len(self.pools)
    def __getitem__(self, key): return self.pools[key]

    @property
    def first(self): return self.pools[0]

    def __init__(self, name, pools, out_dir):
        # Attributes #
        self.name = name
        self.pools = pools
        self.loaded = False
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def load(self):
        # Children #
        for p in self.pools: p.load()
        # Analysis #
        self.analysis = Analysis(self)
        # Save state #
        self.loaded = True

    def run_pools(self, steps=None, **kwargs):
        # Check loaded #
        for p in self.pools:
            if not p.loaded: p.load()
        # Call function #
        for p in self.pools: p()

    def run_pools_slurm(self, steps=None, **kwargs):
        # Check loaded #
        for p in self.pools:
            if not p.loaded: p.load()
        # Call function #
        for p in self.pools: p.run_slurm(steps, **kwargs)

    def make_plots(self):
        # Check loaded #
        if not self.loaded: self.load()
        for cls_name in aggregate_plots.__all__:
            cls = getattr(aggregate_plots, cls_name)
            cls(self).plot()

    def make_zipfile(self):
        # Delete current report #
        if os.path.exists(self.p.report_zip): shutil.rmtree(self.p.report_zip)
        report_dir = self.base_dir + 'report/'
        if os.path.exists(report_dir): shutil.rmtree(report_dir)
        # Copy the experiment results dir #
        shutil.copytree(self.p.results_dir, report_dir)
        # Add the pool results #
        for p in self.pools:
            shutil.copytree(p.p.graphs_dir, report_dir + 'pool%i/' % p.num)
            #for g in l.groups:
            #    shutil.copytree(g.p.graphs_dir, report_dir + 'pool%i/' % l.num + g.short_name + '/')
        # Zip it #
        sh.tar('-zc', '-C', self.base_dir, '-f', self.p.report_zip, 'report')
        shutil.rmtree(report_dir)
