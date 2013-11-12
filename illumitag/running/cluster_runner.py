# Built-in modules #

# Internal modules #
from illumitag.running import Runner
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class ClusterRunner(Runner):
    """Will run stuff on an cluster"""
    default_time = '1-00:00:00'

    default_steps = [
        ### Start ###
        {'combine_reads':             {}},
        ### OTUs ###
        {'run_uparse':                {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.cluster = parent, parent
        self.samples = parent.samples

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        command = """steps = %s
                     names = %s
                     samples = [s for s in illumitag.samples if s.name in names][0]
                     cluster = Cluster(samples)
                     cluster.run(steps)""" % (steps, [s.name for s in self.samples])
        # Test case #
        if self.analysis.aggregate.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        self.slurm_job = SLURMJob(command, self.cluster.p.logs_dir,
                                  job_name="illumitag_cluster_%i" % len(self.samples),
                                  **kwargs)
        self.slurm_job.launch()
