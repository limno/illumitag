b'This module needs Python 2.6 or later.'

# Special variables #
__version__ = '0.1.0'

# Built-in modules #
import os, sys, glob

# Internal modules #
from pools import Pool
from runs import Run

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
# Output directory #
out_dir = home + 'ILLUMITAG/views/pools/'

# Get pool files #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)
project_dir = os.path.abspath(module_dir + '/../') + '/'
pools_dir = project_dir + 'pools/'

# Load all pools #
json_paths = glob.glob(pools_dir + '*.json')
pools = [Pool(j, out_dir) for j in json_paths]

# Compose in to runs #
run_nums = list(set([p.run_num for p in pools])).sort()
runs = [Run(num, [p for p in pools if p.run_num==num]) for num in run_nums]
for p in pools: p.run = [r for r in runs if r.num == p.run_num][0]