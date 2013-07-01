# Built-in modules #

# Internal modules #
from illumitag.common import flatten

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
        pass