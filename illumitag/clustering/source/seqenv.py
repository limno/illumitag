# Built-in modules #

# Internal modules #
from illumitag.common.autopaths import AutoPaths

# Third party modules #
import pandas

###############################################################################
class Seqenv(object):
    """Base class for Seqenv results processing"""

    all_paths = """
    /centers_N5000_blast_F_ENVO_OTUs.csv
    /centers_N5000_blast_F_ENVO_OTUs_labels.csv
    """

    def __init__(self, parent, base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        # Inherited #
        self.samples = self.parent.samples
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.seqenv
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def frame(self):
        return pandas.io.parsers.read_csv(self.p.labels, sep=',', index_col=0, encoding='utf-8')