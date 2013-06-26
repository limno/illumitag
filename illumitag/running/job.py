# Built-in modules #

# Internal modules #
from common import flatten

# Third party modules #

# Constants #

###############################################################################
class Job(object):
    """A bunch of things to do with the data"""

    def __init__(self, pools=None, runs=None, samples=None, steps=None, slurm=False):
        if runs:
            self.pools = flatten([r.pools for r in runs])
        if pools:
            self.pools = pools
