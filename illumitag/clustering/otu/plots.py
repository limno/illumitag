# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import split_thousands

# Third party modules #
from matplotlib import pyplot
import brewer2mpl

# Constants #
__all__ = ['ClusterSizes', 'PresencePerSample', 'PresencePerOTU', 'TaxaBarstack']

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
        fig.suptitle('Clustering method: %s' % self.parent.title)
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
        fig.suptitle('Clustering method: %s' % self.parent.title)
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
        fig.suptitle('Clustering method: %s' % self.parent.title)
        axes.set_xlabel('Number of samples an OTU appears in (max. %i)' % self.parent.otu_table.shape[0])
        axes.set_ylabel('Number of OTUs that appear in these many samples')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class TaxaBarstack(Graph):
    """Distribution of named species by sample"""
    short_name = 'taxa_barstack'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Colors #
        colors = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
        colors.reverse()
        colors += brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
        colors += brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Species relative abundances per sample')
        axes.set_ylabel('Relative abundances in percent')
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=5)
        # Save it #
        self.save_plot(fig, axes, width=24.0, height=14.0, bottom=0.27, top=0.97, left=0.04, right=0.98)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)
