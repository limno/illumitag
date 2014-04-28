# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.graphs.hierarchical_heatmap import HiearchicalHeatmap
from illumitag.common import split_thousands

# Third party modules #
import matplotlib
from matplotlib import pyplot

# Constants #
__all__ = ['ClusterSizes', 'PresencePerSample', 'PresencePerOTU', 'CumulativePresence', 'CumulativePresenceScaled', 'HeatmapOTU']

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
    """Histogram of otu sizes (in terms of number of samples present in them)"""
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

################################################################################
class CumulativePresence(Graph):
    """Cumulative graph of cluster presence in samples. This means something such as:
    - 0% of OTUs appear in 100% of the samples,
    - 10% of OTUs appear in 90% of the samples,
    - 90% of OTUs appear in 1% of the samples"""
    short_name = 'cumulative_presence'

    def plot(self):
        # Number of samples #
        num_of_otus = self.parent.otu_table.shape[1]
        num_of_samples = self.parent.otu_table.shape[0]
        samples_index = list(reversed(range(1,num_of_samples+1)))
        # Get value frequencies #
        counts = self.parent.otu_table.astype(bool).sum(axis=0).value_counts()
        # Add missing values #
        for n in samples_index:
            if n not in counts:
                counts = counts.set_value(n,0)
        # Sort it #
        counts = counts.sort_index(ascending=False)
        # Cumulative sum #
        self.y = list(counts.cumsum())
        # Percentage of samples #
        self.x = [n/num_of_samples for n in samples_index]
        # Make step plot #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.step(self.x, self.y, fillstyle='bottom')
        axes.set_title('Cumulative graph of OTU presence in samples for %s OTUs' % num_of_otus)
        fig.suptitle('Clustering method: %s' % self.parent.otu.title)
        axes.set_xlabel('Fraction of samples (100%% equates %i samples)' % num_of_samples)
        axes.set_ylabel('Number of OTUs that appear in this fraction of samples or more')
        axes.invert_xaxis()
        axes.set_xticks([min(self.x) + (max(self.x)-min(self.x))* n / 9 for n in range(10)])
        axes.set_yscale('log')
        axes.xaxis.grid(True)
        # Set percentage #
        percentage = lambda x, pos: '%1.0f%%' % (x*100.0)
        axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        # Save it #
        self.save_plot(fig, axes)
        pyplot.close(fig)

################################################################################
class CumulativePresenceScaled(Graph):
    """The sample thing but scaled by the size of the OTU"""
    short_name = 'cumulative_presence_scaled'

    def plot(self):
        # Number of samples #
        num_of_samples = self.parent.otu_table.shape[0]
        samples_index = list(reversed(range(1,num_of_samples+1)))
        presences = self.parent.otu_table.astype(bool).sum(axis=0)
        reads_per_otu = self.parent.otu_table.sum(axis=0)
        # Make values #
        self.y = []
        count = 0
        for n in samples_index:
            count += reads_per_otu[presences[presences==n].index].sum()
            self.y.append(count)
        # Percentage of samples #
        self.x = [n/num_of_samples for n in samples_index]
        # Make step plot #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.step(self.x, self.y)
        axes.set_title('Cumulative graph of OTU presence in samples scaled by number of reads')
        fig.suptitle('Clustering method: %s' % self.parent.otu.title)
        axes.set_xlabel('Fraction of samples (100%% equates %i samples)' % num_of_samples)
        axes.set_ylabel('Number of reads in all OTUs that appear in these many samples or more.')
        axes.invert_xaxis()
        axes.xaxis.grid(True)
        # Set percentage #
        percentage = lambda x, pos: '%1.0f%%' % (x*100.0)
        axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        # Save it #
        self.save_plot(fig, axes, left=0.1, sep=('y'))
        pyplot.close(fig)

################################################################################
class HeatmapOTU(HiearchicalHeatmap):
    """Heatmap of the whole OTU table, with hierarchical clustering"""
    short_name = 'otu_heatmap'

    def plot(self):
        # Reference #
        self.frame = self.parent.otu_table.transpose()
        self.row_method = None
        self.fig_height = 200

        # Super #
        fig, axm, axcb, cb = super(HeatmapOTU, self).plot()
        cb.set_label("Number of sequences in cluster")
        # Save #
        pyplot.savefig(self.path)
