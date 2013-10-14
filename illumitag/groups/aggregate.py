# Futures #
from __future__ import division

# Built-in modules #
import os, shutil, re

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.common import slurm
from illumitag.graphs import aggregate_plots
from illumitag.analysis import Analysis
from illumitag.reporting import Reporter
from illumitag.fasta.other import CollectionPairedFASTQ

# Third party modules #
import sh, pandas
from dateutil.parser import parse as dateutil_parse

###############################################################################
class Collection(object):
    """A collection of aggregates."""

    def __repr__(self): return 'Collection: %s' % (self.children)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)

    def __init__(self, children):
        self.children = children

    @property
    def first(self): return self.children[0]

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return [c for c in self.children if c.name == key.lower()][0]
        elif isinstance(key, int):
            if hasattr(self.first, 'num'): return [c for c in self.children if c.num == key][0]
            else: return self.children[key]
        else: raise TypeError('key')

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several pools."""

    all_paths = """
    /graphs/
    /logs/
    /results/slurm_report.csv
    """

    def __repr__(self): return '<%s object "%s" with %i pools>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.pools)
    def __len__(self): return len(self.pools)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int): return self.children[key]
        else: raise TypeError('key')

    @property
    def first(self): return self.pools[0]

    @property
    def count(self): return sum(map(len, self.outcomes))

    def __init__(self, name, pools, out_dir):
        # Attributes #
        self.name = name
        self.pools = pools
        self.loaded = False
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def load(self):
        # Children #
        for p in self.pools:
            if not p.loaded: p.load()
        # Make Outcomes #
        self.no_barcodes   = CollectionPairedFASTQ(p.no_barcodes for p in self)
        self.one_barcodes  = CollectionPairedFASTQ(p.one_barcodes for p in self)
        self.same_barcodes = CollectionPairedFASTQ(p.same_barcodes for p in self)
        self.bad_barcodes  = CollectionPairedFASTQ(p.bad_barcodes for p in self)
        self.good_barcodes = CollectionPairedFASTQ(p.good_barcodes for p in self)
        self.outcomes = (self.good_barcodes, self.no_barcodes, self.one_barcodes, self.same_barcodes, self.bad_barcodes)
        # Analysis #
        self.analysis = Analysis(self)
        # Reporting #
        self.reporter = Reporter(self)
        # Graphs #
        self.graphs = [getattr(aggregate_plots, cls_name)(self) for cls_name in aggregate_plots.__all__]
        # Save state #
        self.loaded = True

    def run_pools(self, steps=None, **kwargs):
        # Check loaded #
        for p in self.pools:
            if not p.loaded: p.load()
        # Call function #
        for p in self.pools: p()

    def run_pools_slurm(self, steps=None, **kwargs):
        return [p.run_slurm(steps, **kwargs) for p in self.pools]

    def make_plots(self):
        if not self.loaded: self.load()
        for graph in self.graphs: graph.plot()

    def make_slurm_report(self):
        running_jobs_names = [j['name'] for j in slurm.running_jobs_info()]
        queued_jobs_names = [j['name'] for j in slurm.running_jobs_info()]
        for p in self:
            # Loaded #
            if not p.loaded: p.load()
            # Let's see if it didn't fail #
            p.job_state = "Failed"
            # Running #
            if str(p) in queued_jobs_names:
                p.job_state = "Queued"
                continue
            if str(p) in running_jobs_names:
                p.job_state = "Running"
                continue
            # Out path #
            job_out_path = p.runner.latest_log + 'run.out'
            # No file #
            if not os.path.exists(job_out_path):
                p.job_state = "Missing"
                continue
            # Read the log #
            with open(job_out_path, 'r') as handle: job_out = handle.read()
            # Success #
            if 'Success' in job_out: p.job_state = "Success"
            # Problems #
            if 'CANCELLED' in job_out: p.job_state = "Slurm CANCELLED"
            if '(core dumped)' in job_out: p.job_state = "Core DUMPED"
            # Start and end time #
            start_time = re.findall('^SLURM: start at (.+) on .+$', job_out, re.M)
            if start_time: p.runner.job_start_time = dateutil_parse(start_time[0])
            end_time = re.findall('^SLURM: end at (.+)$', job_out, re.M)
            if end_time: p.runner.job_end_time = dateutil_parse(end_time[0])
            # Total time #
            if start_time and end_time:
                p.runner.job_runtime = p.runner.job_end_time - p.runner.job_start_time
        # Make report #
        rows = [str(p) for p in self]
        columns = ['Name', 'Project', 'State', 'Start time', 'End time', 'Run time']
        data = [(p.long_name, p.project.name, p.job_state,
                str(p.runner.job_start_time), str(p.runner.job_end_time), str(p.runner.job_runtime))
                for p in self]
        frame = pandas.DataFrame(data, index=rows, columns=columns)
        frame.to_csv(self.p.slurm_report)

    def make_zipfile(self):
        # Delete current report #
        if os.path.exists(self.p.report_zip): shutil.rmtree(self.p.report_zip)
        report_dir = self.base_dir + 'report/'
        if os.path.exists(report_dir): shutil.rmtree(report_dir)
        # Copy the experiment results dir #
        shutil.copytree(self.p.results_dir, report_dir)
        # Add the pool results #
        for p in self.pools:
            shutil.copytree(p.p.graphs_dir, report_dir + 'pool%i/' % p.num)
            #for g in l.groups:
            #    shutil.copytree(g.p.graphs_dir, report_dir + 'pool%i/' % l.num + g.short_name + '/')
        # Zip it #
        sh.tar('-zc', '-C', self.base_dir, '-f', self.p.report_zip, 'report')
        shutil.rmtree(report_dir)
