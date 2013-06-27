# Built-in modules #

# Internal modules #
from illumitag.common import flatten
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
def make_jobs(pools=None, runs=None, steps=None):
    if runs:
        pools = flatten([r.pools for r in runs])
    if pools:
        pools = pools
    return [PoolJob(p, steps) for p in pools]

###############################################################################
class PoolJob(object):
    def __init__(self, pool, steps=None):
        self.pool = pool
        self.steps = steps

    def run(self, steps, **kwargs):
        # Make script #
        command = """
            steps = %s
            pool = [p for p in illumitag.pools if str(p)==%s][0]
            pool(steps)
        """
        command = command % (steps, self.pool)
        command = '\n'.join(l.lstrip(' ') for l in command.split('\n') if l)
        # Send it #
        slurm_job = SLURMJob(command, log_dir=self.pool.p.logs_dir, **kwargs)
        slurm_job.launch()