# Built-in modules #

# Internal modules #
import illumitag
from illumitag.groups.aggregate import Collection, Aggregate
from illumitag.common import AutoPaths
from illumitag.common.slurm import SLURMJob

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

    def run_analysis_slurm(self, steps=None, **kwargs):
        # Check loaded #
        if not self.loaded: self.load()
        # Make script #
        command = """steps = %s
                     proj = [pj for pj in illumitag.projects if pj.name=='%s'][0]
                     proj.load()
                     proj.analysis.run()""" % (steps, self.name)
        # Test case #
        if self.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = '1-00:00:00'
        if 'email' not in kwargs: kwargs['email'] = None
        self.slurm_job = SLURMJob(command, self.p.logs_dir, job_name="illumitag_"+self.name, **kwargs)
        self.slurm_job.launch()