# Built-in modules #

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.clustering.statistics.nmds import NMDS
from illumitag.clustering.statistics.permanova import PERMANOVA
from illumitag.clustering.statistics.betadis import BetaDispersion

# Third party modules #

###############################################################################
class StatsOnOTUs(object):

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    """

    def __init__(self, parent):
        # Save parent #
        self.otu, self.parent = parent, parent
        # Paths #
        self.p = AutoPaths(self.parent.p.stats_dir, self.all_paths)
        # Children #
        self.nmds = NMDS(self, self.parent.otu_csv)
        self.permanova = PERMANOVA(self)
        self.betadis = BetaDispersion(self)

    def run(self):
        self.nmds.run()
        self.permanova.run()
        self.betadis.run()

###############################################################################
class StatsOnComposition(object):

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    """

    def __init__(self, parent):
        # Save parent #
        self.composition, self.parent = parent, parent
        # Paths #
        self.p = AutoPaths(self.parent.p.stats_dir, self.all_paths)
        # Children #
        self.nmds = NMDS(self, self.parent.taxa_csv)

    def run(self):
        self.nmds.run()
