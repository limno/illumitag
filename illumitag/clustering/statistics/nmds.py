# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.common.conversion import r_matrix_to_dataframe, pandas_df_to_r_df
from illumitag.graphs import Graph

# Third party modules #
from rpy2 import robjects as ro
from matplotlib import pyplot

################################################################################
class GraphNMDS(Graph):
    """Non-metric dimensional scaling"""
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
            pyplot.annotate(names[i], size=9, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes, bottom=0.03, top=0.97)
        pyplot.close(fig)

###############################################################################
class NMDS(object):

    all_paths = """
    /lorem
    """

    def __init__(self, parent, csv, calc_distance=True):
        # Save parent #
        self.stat, self.parent = parent, parent
        self.csv = csv
        # Options #
        self.calc_distance = calc_distance
        # Paths #
        self.base_dir = self.parent.p.nmds_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Graph #
        self.graph = GraphNMDS(self, base_dir=self.base_dir)

    def run(self):
        # Load dataframe #
        ro.r("library(vegan)")
        ro.r("table = read.table('%s', sep='\t', header=TRUE, row.names='X')" % (self.csv))
        # Run computation #
        if self.calc_distance: ro.r("nmds = metaMDS(table, distance='horn', trymax=200)")
        else:                  ro.r("nmds = metaMDS(table, trymax=200)")
        # Extract result #
        ro.r("coord = scores(nmds)")
        ro.r("loadings = nmds$species")
        # Retrieve values #
        self.coords = r_matrix_to_dataframe(ro.r.coord)
        # No loadings without distance #
        if self.calc_distance: self.loadings = r_matrix_to_dataframe(ro.r.loadings)
        else:                  self.loadings = False
        # Plot it #
        self.graph.plot()

    def run_df(self):
        """Unfortunately this doesn't seem to work (yet)"""
        # Get frame #
        self.frame = self.parent.parent.frame
        # Call R #
        rframe = pandas_df_to_r_df(self.frame)
        ro.r("library(vegan)")
        nmds = ro.r['metaMDS'](rframe, distance='horn', trymax=200)
        # Retrieve values #
        self.coords = r_matrix_to_dataframe(ro.r['scores'](nmds))
        self.loadings = list(nmds.rx2('species'))
        # Plot it #
        self.graph.plot()