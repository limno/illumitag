b"""This module needs Python 2.6 or later."""

# Special variables #
__version__ = '0.1.0'

# Built-in modules #
import os, sys, glob

# Internal modules #
from illumitag.groups.pools import Pool
from illumitag.groups.runs import Run
from illumitag.groups.projects import Project
from illumitag.common import dependencies

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
# Check dependencies #
dependencies.check_modules()
dependencies.check_executables()

# Output directory #
view_dir = out_dir = home + 'ILLUMITAG/views/'

# Get pool files #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)
project_dir = os.path.abspath(module_dir + '/../') + '/'
pools_dir = project_dir + 'pools/'

# Load all pools #
json_paths = glob.glob(pools_dir + '*.json')
pools = [Pool(j, view_dir + 'pools/') for j in json_paths]

# Compose into runs #
run_nums = sorted(list(set([p.run_num for p in pools])))
runs = [Run(num, [p for p in pools if p.run_num==num], view_dir + 'runs/') for num in run_nums]
for p in pools: p.run = [r for r in runs if r.num == p.run_num][0]

# Compose into projects #
proj_names = sorted(list(set([p.project_short_name for p in pools])))
projects = [Project(name, [p for p in pools if p.project_short_name==name], view_dir + 'projects/') for name in proj_names]
for pool in pools: p.project = [proj for proj in projects if proj.name == pool.project_short_name][0]
