# Built-in modules #

# Internal modules #

# Third party modules #

###############################################################################
class OTUs(object):
    """Base class for OTU making"""

    def make_plots(self):
        for graph in self.graphs: graph.plot()