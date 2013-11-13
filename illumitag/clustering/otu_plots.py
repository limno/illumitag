# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import split_thousands

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['ClusterDistribution']

################################################################################
class ClusterDistribution(Graph):
    """Distribution of cluster sizes"""
    short_name = 'cluster_distribution'

    def plot(self):
        # Column sums and count frequencies #
        distrib = self.parent.frame.sum().value_counts()
        x = distrib.keys()
        y = distrib.values
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_xscale('symlog')
        axes.set_yscale('log')
        axes.set_title('Distribution of cluster sizes for %s clusters' % split_thousands(sum(y)))
        fig.suptitle('Clustering method: %s' % self.parent.title)
        axes.set_xlabel('Total children of a cluster')
        axes.set_ylabel('Number of clusters with this many children')
        axes.xaxis.grid(False)
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)