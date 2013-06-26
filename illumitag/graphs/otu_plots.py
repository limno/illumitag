# Futures #
from __future__ import division

# Built-in modules #
import time, getpass, csv
from collections import defaultdict, OrderedDict
from itertools import izip

# Internal modules #
from illumitaq.util import save_plot
from illumitaq.common import split_thousands

# Third party modules #
import pandas
from matplotlib import pyplot

###############################################################################
def cluster_distribution(otus):
    """Distribution of cluster sizes in all samples"""
    # Read data #
    dist = defaultdict(int)
    handle = open(otus.p.clusters_otus_txt)
    for line in handle: dist[len(line.split())-1] += 1
    dist = OrderedDict(sorted(dist.items()))
    # Make scatter #
    fig = pyplot.figure()
    axes = fig.add_subplot(111)
    x = dist.keys()
    y = dist.values()
    axes.plot(x, y, 'ro')
    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_title('Distribution of cluster sizes for %s clusters' % split_thousands(sum(y)))
    fig.suptitle('Clustering method: %s' % otus.method)
    axes.set_xlabel('Total children of a cluster')
    axes.set_ylabel('Frequency of cluster size')
    axes.xaxis.grid(False)
    # Add annotations #
    for i in range(5):
        pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                        textcoords = 'offset points', ha = 'left', va = 'center',
                        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
    # Save it #
    save_plot(fig, axes, otus.p.cluster_hist_pdf)

###############################################################################
def otu_distribution(otus):
    """Distribution of OTUS counts in all samples"""
    # Read data #
    dist = defaultdict(int)
    handle = open(otus.p.csv_table)
    handle.next()
    for line in handle: dist[sum([1 for x in line.split()[1:] if x != '0.0'])] += 1
    dist = OrderedDict(sorted(dist.items()))
    # Make scatter #
    fig = pyplot.figure()
    axes = fig.add_subplot(111)
    x = dist.keys()
    y = dist.values()
    axes.plot(x, y, 'ro')
    axes.set_yscale('log')
    axes.set_title('Distribution of OTU presence/absence counts for %s OTUs' % split_thousands(sum(y)))
    fig.suptitle('Clustering method: %s' % otus.method)
    axes.set_xlabel('Number of samples a given OTU appears in')
    axes.set_ylabel('Frequency of OTU span')
    axes.xaxis.grid(False)
    # Add annotations #
    for i in range(5):
        pyplot.annotate("%i: %s" % (x[i], split_thousands(y[i])), size=13, xy = (x[i], y[i]), xytext = (10, 0),
                        textcoords = 'offset points', ha = 'left', va = 'center',
                        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
    # Save it #
    save_plot(fig, axes, otus.p.otu_hist_pdf)

###############################################################################
def otu_sums(otus):
    """Distribution of OTUS sums per sample"""
    # Read data #
    table = pandas.read_csv(otus.p.table_trimmed, sep = '\t', index_col=0)
    sums = table.sum(axis=1)
    # Make hist #
    fig = pyplot.figure()
    axes = sums.hist(color='gray', bins=40)
    fig = pyplot.gcf()
    axes.set_title('Distribution of OTU appearance sums per sample')
    fig.suptitle('Clustering method: %s' % otus.method)
    axes.set_xlabel('Sum of all OTU appearances in a given sample')
    axes.set_ylabel('Frequency of sample size')
    axes.xaxis.grid(False)
    # Save it #
    save_plot(fig, axes, otus.p.otu_sums_hist_pdf, sep=('x'))

def sample_sums(otus):
    """Distribution of OTUS sums per sample"""
    # Read data #
    table = pandas.read_csv(otus.p.table_trimmed, sep = '\t', index_col=0)
    sums = table.sum(axis=0)
    # Make hist #
    fig = pyplot.figure()
    axes = fig.add_subplot(111)
    axes.hist(sums, 40, color='gray', log=True)
    axes.set_title('Distribution of OTU appearance sums per OTU')
    fig.suptitle('Clustering method: %s' % otus.method)
    axes.set_xlabel('Sum of OTU appearance in all samples for a given OTU')
    axes.set_ylabel('Frequency of OTU appearance total')
    axes.xaxis.grid(False)
    # Save it #
    save_plot(fig, axes, otus.p.sample_sums_hist_pdf, sep=('x'))

###############################################################################
def taxa_hist(otus, pools=["pool1", "pool2", "pool3", "pool4", "pool5"]):
    """Distribution of barcodes within a pool"""
    # Make unfiltered version #
    handle = open(otus.p.unfiltered_table, "w")
    rows = izip(*csv.reader(open(otus.p.csv_table), delimiter='\t'))
    csv.writer(handle, delimiter='\t').writerows(rows)
    handle.close()
    # Parse taxonomy #
    otu_to_phylums = {}
    for line in open(otus.p.rep_set_taxonomy):
        otu, read, taxa = line.split()
        taxa = taxa.split(';')
        otu_to_phylums[otu] = taxa[1]
    # For every pool #
    for pool in pools:
        # Axes #
        columns = set(otu_to_phylums.values())
        rows = [pool + '-' + name for name in otus.procedure.exp1.first.barcodes.names]
        # Parse OTU table #
        handle = open(otus.p.unfiltered_table)
        otu_names = handle.next().split()[1:]
        # Main loop #
        sample_to_breakdown = {}
        for line in handle:
            if line[0:5] != pool: continue
            line = line.split()
            sample, otu_counts = line[0], line[1:]
            breakdown = defaultdict(int)
            for otu, count in izip(otu_names, otu_counts):
                if count == '0.0': continue
                breakdown[otu_to_phylums[otu]] += int(float(count))
            total = sum(breakdown.values())
            for key in breakdown: breakdown[key] = (breakdown[key] / total) * 100.0
            sample_to_breakdown[sample] = breakdown
        # Frame #
        data = [[sample_to_breakdown[sample][phyla] for phyla in columns] for sample in rows]
        otus.taxa_hist = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = otus.taxa_hist.plot(kind='bar', stacked=True)
        fig = pyplot.gcf()
        axes.set_title('Distribution of taxa by sample for %s' % pool)
        axes.set_ylabel('Fraction of reads classified at phylum level')
        axes.xaxis.grid(False)
        # Size #
        fig.set_figwidth(18.0)
        fig.set_figheight(10.0)
        #fig.subplots_adjust(hspace=0.0, bottom=0.1, top=0.93, left=0.06, right=0.98)
        # Data and source #
        fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
        fig.text(0.01, 0.98, 'user: ' + getpass.getuser(), horizontalalignment='left')
        # Legend #
        box = axes.get_position()
        axes.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Save it #
        fig.savefig(getattr(otus.p, pool + '_pdf'))
        otus.taxa_hist.to_csv(getattr(otus.p, pool + '_csv'))
