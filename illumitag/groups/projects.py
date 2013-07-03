# Built-in modules #

# Internal modules #
import illumitag
from illumitag.groups.aggregate import Collection, Aggregate
from illumitag.common import AutoPaths

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several pools possibly spanning several runs."""

    def __repr__(self): return '<%s object "%s" with %i pools>' % \
                               (self.__class__.__name__, self.name, len(self))

    @property
    def long_name(self): return self.first.project_long_name

    def __init__(self, name, pools, projs_dir):
        # Attributes #
        self.name = name
        self.pools = pools
        self.loaded = False
        # Dir #
        self.base_dir = projs_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra #
        self.meta_data_path = illumitag.repos_dir + 'projects/' + self.name + '.csv'

    def run_pools_slurm(self, steps=None, **kwargs):
        # Test case #
        if self.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Call function #
        return Aggregate.run_pools_slurm(self, **kwargs)

    def run_analysis_slurm(self, *args, **kwargs):
        if not self.loaded: self.load()
        self.analysis.run_slurm(*args, **kwargs)