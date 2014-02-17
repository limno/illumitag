b"""This module needs Python 2.7.x"""

# Special variables #
__version__ = '0.2.0'

# Built-in modules #
import os, sys, glob

# Internal modules #
from illumitag.groups.pools import Pool
from illumitag.groups.runs import Runs, Run
from illumitag.groups.aggregate import Aggregate
from illumitag.groups.projects import Projects, Project
from illumitag.groups.presamples import Presample
from illumitag.groups.pyrosample import Pyrosample
from illumitag.common import dependencies

# Constants #
home = os.environ['HOME'] + '/'

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

###############################################################################
# Check dependencies #
dependencies.check_modules()
dependencies.check_executables()

# Output directory #
view_dir = out_dir = home + 'ILLUMITAG/views/'

# Get pool files #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)
repos_dir = os.path.abspath(module_dir + '/../') + '/'
pools_dir = repos_dir + 'json/pools/'

# Load all pools #
json_paths = glob.glob(pools_dir + '*.json')
pools = [Pool(j, view_dir + 'pools/') for j in json_paths]
pools.sort(key=lambda x: str(x))

# Compose into runs #
run_nums = sorted(list(set([p.run_num for p in pools]))) # [1,2,3,4,5]
runs = [Run(num, [p for p in pools if p.run_num==num], view_dir + 'runs/') for num in run_nums]
runs = Runs(runs)
for p in pools: p.run = runs[p.run_num]

# Compose into projects #
proj_names = sorted(list(set([p.project_short_name for p in pools])))
projects = [Project(name, [p for p in pools if p.project_short_name==name], view_dir + 'projects/') for name in proj_names]
projects = Projects(projects)
for p in pools: p.project = projects[p.project_short_name]

# Make an aggregate with all pools #
aggregate = Aggregate('all', pools, view_dir + 'aggregates/')

###############################################################################
# Load all presamples #
presamples_dir = repos_dir + 'json/presamples/'
json_paths = glob.glob(presamples_dir + '*.json')
presamples = [Presample(j, view_dir + 'presamples/') for j in json_paths]
presamples.sort(key=lambda x: str(x))

# Load all legacy pyrosamples #
presamples_dir = repos_dir + 'json/pyrosample/'
json_paths = glob.glob(presamples_dir + '*.json')
presamples = [Pyrosample(j, view_dir + 'pyrosamples/') for j in json_paths]
presamples.sort(key=lambda x: str(x))
