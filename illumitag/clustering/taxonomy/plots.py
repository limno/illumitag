# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import split_thousands

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['ClusterSizes', 'PresencePerSample', 'PresencePerOTU']

################################################################################
class ClusterSizes(Graph):
    """Distribution of cluster sizes in loglog"""
    short_name = 'cluster_sizes'

    def plot(self):
        # Sum by column and count frequencies #
        distrib = self.parent.otu_table.sum().value_counts()
        x = distrib.keys()
        y = distrib.values
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.plot(x, y, 'ro')
        axes.set_xscale('symlog')
        axes.set_yscale('log')
        axes.set_title('Distribution of cluster sizes for %s clusters' % split_thousands(sum(y)))
        fig.suptitle('Clustering method: %s' % self.parent.otu.title)
        axes.set_xlabel('Number of children in a cluster')
        axes.set_ylabel('Number of clusters with this many children')
        axes.xaxis.grid(False)
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)
        pyplot.close(fig)

################################################################################
class PresencePerSample(Graph):
    """Histogram of sample sizes (in terms of number of OTUs present in them)"""
    short_name = 'sample_sums'

    def plot(self):
        # Sum by row and count frequencies #
        self.frame = self.parent.otu_table.astype(bool).sum(axis=1)
        # Make hist #
        fig = pyplot.figure()
        axes = self.frame.hist(color='gray', bins=40)
        fig = pyplot.gcf()
        axes.set_title('Histogram of OTU appearance sums per sample')
        fig.suptitle('Clustering method: %s' % self.parent.otu.title)
        axes.set_xlabel('Number of OTUs present (non-null) in a sample')
        axes.set_ylabel('Number of samples with this many OTUs in them')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class PresencePerOTU(Graph):
    """Histogram of sample sizes (in terms of number of OTUs present in them)"""
    short_name = 'otu_sums'

    def plot(self):
        # Sum by row and count frequencies #
        self.frame = self.parent.otu_table.astype(bool).sum(axis=0)
        # Make hist #
        fig = pyplot.figure()
        axes = self.frame.hist(color='gray', bins=40)
        fig = pyplot.gcf()
        axes.set_title('Histogram of OTU appearance sums per OTU')
        fig.suptitle('Clustering method: %s' % self.parent.otu.title)
        axes.set_xlabel('Number of samples an OTU appears in (max. %i)' % self.parent.otu_table.shape[0])
        axes.set_ylabel('Number of OTUs that appear in these many samples')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)