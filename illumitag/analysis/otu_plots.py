# Futures #
from __future__ import division
from collections import defaultdict, OrderedDict

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import split_thousands

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['ClusterDistributionRaw', 'ClusterDistributionFiltered', 'OTUDistribution', 'OTUSums', 'SampleSums']


################################################################################
class ClusterDistributionRaw(Graph):
    """Distribution of cluster sizes within a pool after clustering"""
    short_name = 'cluster_distribution_raw'

    def plot(self):
        dist = defaultdict(int)
        handle = open(self.parent.p.clusters_otus_txt)
        for line in handle: dist[len(line.split())-1] += 1
        dist = OrderedDict(sorted(dist.items()))
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        x = dist.keys()
        y = dist.values()
        axes.plot(x, y, 'ro')
        axes.set_xscale('symlog')
        axes.set_yscale('log')
        axes.set_title('Distribution of cluster sizes for %s clusters' % split_thousands(sum(y)))
        fig.suptitle('Clustering method: %s' % self.parent.method)
        axes.set_xlabel('Total children of a cluster')
        axes.set_ylabel('Frequency of cluster size')
        axes.xaxis.grid(False)
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)

################################################################################
class ClusterDistributionFiltered(Graph):
    """Distribution of barcodes within a pool"""
    short_name = 'cluster_distribution_filtered'

    def plot(self):
        # Read data #
        dist = defaultdict(int)
        handle = open(self.parent.table_transposed.path)
        handle.next()
        for line in handle: dist[sum([int(x) for x in line.split()[1:] if x != '0'])] += 1
        # Zero is Nan on log plot #
        #if 0 in dist: dist[0.1] = dist.pop(0)
        dist = OrderedDict(sorted(dist.items()))
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        x = dist.keys()
        y = dist.values()
        axes.plot(x, y, 'ro')
        axes.set_xscale('symlog')
        axes.set_yscale('log')
        axes.set_title('Distribution of cluster sizes for %s clusters (after filtering)' % split_thousands(sum(y)))
        fig.suptitle('Clustering method: %s' % self.parent.method)
        axes.set_xlabel('Total children of a cluster')
        axes.set_ylabel('Frequency of cluster size')
        axes.xaxis.grid(False)
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)

###############################################################################
class OTUDistribution(Graph):
    """Distribution of OTUS counts in all samples"""
    short_name = 'otu_distribution'

    def plot(self):
        # Read data #
        dist = defaultdict(int)
        handle = open(self.parent.table.path)
        handle.next()
        for line in handle: dist[sum([1 for x in line.split()[1:] if x != '0'])] += 1
        dist = OrderedDict(sorted(dist.items()))
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        x = dist.keys()
        y = dist.values()
        axes.plot(x, y, 'ro')
        axes.set_yscale('log')
        axes.set_title('Distribution of OTU presence/absence counts for %s OTUs' % split_thousands(sum(y)))
        fig.suptitle('Clustering method: %s' % self.parent.method)
        axes.set_xlabel('Number of samples a given OTU appears in')
        axes.set_ylabel('Frequency of OTU span')
        axes.xaxis.grid(False)
        # Add annotations #
        for i in range(min(5,len(x))):
            pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                            textcoords = 'offset points', ha = 'left', va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)

###############################################################################
class OTUSums(Graph):
    """Distribution of OTUS sums per sample"""
    short_name = 'otu_sums'

    def plot(self):
        # Read data #
        table = pandas.read_csv(self.parent.table_filtered.path, sep = '\t', index_col=0)
        sums = table.sum(axis=1)
        # Make hist #
        fig = pyplot.figure()
        axes = sums.hist(color='gray', bins=40)
        fig = pyplot.gcf()
        axes.set_title('Distribution of OTU appearance sums per sample')
        fig.suptitle('Clustering method: %s' % self.parent.method)
        axes.set_xlabel('Sum of all OTU appearances in a given sample')
        axes.set_ylabel('Frequency of sample size')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))

###############################################################################
class SampleSums(Graph):
    """Distribution of OTUS sums per OTU"""
    short_name = 'sample_sums'

    def plot(self):
        # Read data #
        table = pandas.read_csv(self.parent.table_filtered.path, sep = '\t', index_col=0)
        sums = table.sum(axis=0)
        # Make hist #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.hist(sums, 40, color='gray', log=True)
        axes.set_title('Distribution of OTU appearance sums per OTU')
        fig.suptitle('Clustering method: %s' % self.parent.method)
        axes.set_xlabel('Sum of OTU appearance in all samples for a given OTU')
        axes.set_ylabel('Frequency of OTU appearance total')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
