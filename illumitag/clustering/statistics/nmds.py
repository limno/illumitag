# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.common.conversion import r_matrix_to_dataframe
from illumitag.graphs import Graph

# Third party modules #
from rpy2 import robjects as ro
from matplotlib import pyplot

################################################################################
class GraphNMDS(Graph):
    """Lorem"""
    short_name = 'nmds'

    def plot(self):
        # Coord #
        x = self.parent.coords['NMDS1'].values
        y = self.parent.coords['NMDS2'].values
        names = self.parent.coords['NMDS1'].keys()
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_title('Non-Metric Multidimensional scaling')
        axes.set_xlabel('Dimension 1')
        axes.set_ylabel('Dimension 2')
        # Add annotations #
        for i in range(len(names)):
            pyplot.annotate(names[i], size=6, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes, width=20.0, height=20.0, bottom=0.03, top=0.97)
        pyplot.close(fig)

###############################################################################
class NMDS(object):

    all_paths = """
    /lorem
    """

    def __init__(self, parent):
        # Save parent #
        self.stat, self.parent = parent, parent
        self.otu = parent.otu
        # Paths #
        self.base_dir = self.parent.p.nmds_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graph #
        self.graph = GraphNMDS(self, base_dir=self.base_dir)

    def run(self):
        # Call R #
        ro.r("library(vegan)")
        ro.r("otu_table = read.table('%s', sep='\t', header=TRUE, row.names='X')" % (self.otu.otu_csv))
        ro.r("nmds = metaMDS(otu_table, distance='horn', trymax=200)")
        ro.r("coord = scores(nmds)")
        ro.r("loadings = nmds$species")
        # Retrieve values #
        self.coords = r_matrix_to_dataframe(ro.r.coord)
        self.loadings = r_matrix_to_dataframe(ro.r.loadings)

