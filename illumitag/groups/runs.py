# Built-in modules #
import os
import xml.etree.ElementTree as etree

# Internal modules #
from aggregate import Collection, Aggregate
from illumitag.common.autopaths import AutoPaths

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Runs(Collection):
    """A collection of runs."""
    pass

###############################################################################
class Run(Aggregate):
    """An illumina run containing several pools."""

    def __repr__(self): return '<%s object number %i with %i pools>' % \
                               (self.__class__.__name__, self.num, len(self))

    @property
    def label(self): return self.first.run_label

    def __init__(self, num, pools, out_dir):
        # Attributes #
        self.num = num
        self.name = "run%i" % num
        self.pools, self.children = pools, pools
        self.loaded = False
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra #
        self.meta_data_path = None
        # Illumina report #
        self.xml_report_path = home + "ILLUMITAG/INBOX/%s/report.xml" % self.label

    def parse_report_xml(self):
        tree = etree.parse(self.xml_report_path)
        root = tree.getroot()
        # The stats of the run #
        lane = root.iter('Lane').next()
        self.report_stats = {}
        self.report_stats['fwd'], self.report_stats['rev'] = map(lambda x: dict(x.items()), lane.iter('Read'))
        # For every sample #
        for sample in root.iter('Sample'):
            label = sample.items()[0][1]
            pool = [p for p in self if p.label == label]
            # Common name mix up #
            if not pool: pool = [p for p in self if p.label.replace("Sample_","") == label]
            pool = pool[0]
            # Get the stats #
            pool.report_stats = {}
            pool.report_stats['fwd'], pool.report_stats['rev'] = map(lambda x: dict(x.items()), sample.iter('Read'))