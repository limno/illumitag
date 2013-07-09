# Built-in modules #

# Internal modules #
from illumitag.running import Runner
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class AnalysisRunner(Runner):
    """Will run stuff on an aggregate's analysis"""
    default_time = '1-00:00:00'

    default_steps = [
        ### Start ###
        {'combine_reads':             {}},
        ### OTUs ###
        {'run_denovo':                {}},
        {'run_subsample':             {}},
        #{'run_open_ref':              {}},
        #{'run_progressive':           {}},
        ### Plots ###
        #{'make_pool_plots':           {'threads':False}},
        #{'make_outcome_plots':        {'threads':False}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.analysis = parent, parent
        self.pools = parent.pools

    def run_slurm(self, steps=None, **kwargs):
        # Check project #
        from illumitag.groups.projects import Project
        if not isinstance(self.analysis.aggregate, Project):
            raise Exception("Can only analyze projects via SLURM.")
        # Make script #
        command = """steps = %s
                     proj = [pj for pj in illumitag.projects if pj.name=='%s'][0]
                     proj.load()
                     proj.analysis.run()""" % (steps, self.analysis.aggregate.name)
        # Test case #
        if self.analysis.aggregate.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Dependencies #
        if all([hasattr(p, 'slurm_job') for p in self.pools]):
            deps = 'afterok:' + ':'.join(str(p.slurm_job.id) for p in self.pools)
        else:
            deps = 'afterok:1'
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        self.slurm_job = SLURMJob(command, self.analysis.p.logs_dir,
                                  job_name="illumitag_" + self.analysis.aggregate.name,
                                  dependency=deps, **kwargs)
        self.slurm_job.launch()
