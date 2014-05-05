# Futures #
from __future__ import division

# Built-in modules #
from collections import Counter, OrderedDict

# Internal modules #
import illumitag
from illumitag.graphs import Graph, cool_colors
from illumitag.helper.silvamod import amplified
from illumitag.common import split_thousands

# Third party modules #
import pandas, matplotlib
from matplotlib import pyplot

################################################################################
class BarcodeStack(Graph):
    short_name = 'barcode_stack'
    bottom = 0.1
    top = 0.95
    left = 0.15
    right = 0.95
    formats = ('pdf', 'eps')

    def plot(self):
        # Data #
        rows  = ['Updated chemistry']
        rows += ['Single-step PCR']
        rows += ['Two-step PCR 3']
        rows += ['Two-step PCR 2']
        rows += ['Two-step PCR 1']
        columns = [o.doc for o in self.parent.first.outcomes]
        data = [[o.count for o in p.outcomes] for p in reversed(self.parent.pools)]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh', stacked=True, color=['g','k','y','orange','r'])
        fig = pyplot.gcf()
        axes.set_xlabel('Number of paired reads')
        axes.xaxis.grid(True)
        axes.yaxis.grid(False)
        axes.legend(prop={'size':10})
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class LengthDistribution(Graph):
    short_name = 'length_distribution'
    bottom = 0.1
    top = 0.95
    left = 0.1
    right = 0.9
    formats = ('pdf', 'eps')

    def plot(self):
        # Data #
        counts = sum((p.quality_reads.only_used.lengths for p in self.parent), Counter())
        # Plot #
        fig = pyplot.figure()
        pyplot.bar(counts.keys(), counts.values(), 1.0, color='gray', align='center', label='Reads from soil sample')
        axes = pyplot.gca()
        axes.set_xlabel('Length of sequence in nucleotides')
        axes.set_ylabel('Number of sequences with this length')
        axes.yaxis.grid(True)
        # Add Silvamod lengths on second scale #
        silva_counts = amplified.lengths
        silvas_axes = axes.twinx()
        silvas_axes.plot(silva_counts.keys(), silva_counts.values(), 'r-', label='Sequences from the silvamod database')
        silvas_axes.set_ylabel('Number of sequences from the silvamod database', color='r')
        for tick in silvas_axes.get_yticklabels(): tick.set_color('r')
        # Legends #
        axes.legend(loc='upper left')
        silvas_axes.legend(loc='upper right')
        # Add separator #
        seperate = lambda y,pos: split_thousands(y)
        silvas_axes.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
        # Change ticks #
        import matplotlib.ticker as mticker
        myLocator = mticker.MultipleLocator(10)
        axes.xaxis.set_major_locator(myLocator)
        axes.set_xlim(400, 500)
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        # Save CSV #
        self.frame = pandas.Series(counts.get(i,0) for i in range(max(counts.keys())+1))
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class FractionTaxaBarStack(Graph):
    short_name = 'fraction_taxa_barstack'
    bottom = 0.4
    top = 0.95
    left = 0.1
    right = 0.95
    formats = ('pdf', 'eps')

    def plot(self):
        # Make Frame #
        self.frame = OrderedDict((('%s - %s' % (p,f), getattr(p.fractions, f).rdp.phyla)
                     for f in ('low', 'med', 'big') for p in self.parent.pools))
        self.frame = pandas.DataFrame(self.frame)
        self.frame = self.frame.fillna(0)
        # Rename #
        new_names = {
            u"run001-pool01 - low": "2-step PCR low",
            u"run001-pool02 - low": "2-step PCR low",
            u"run001-pool03 - low": "2-step PCR low",
            u"run001-pool04 - low": "1-step PCR low",
            u"run002-pool01 - low": "New chem low",
            u"run001-pool01 - med": "2-step PCR med",
            u"run001-pool02 - med": "2-step PCR med",
            u"run001-pool03 - med": "2-step PCR med",
            u"run001-pool04 - med": "1-step PCR med",
            u"run002-pool01 - med": "New chem med",
            u"run001-pool01 - big": "2-step PCR high",
            u"run001-pool02 - big": "2-step PCR high",
            u"run001-pool03 - big": "2-step PCR high",
            u"run001-pool04 - big": "1-step PCR high",
            u"run002-pool01 - big": "New chem high",
        }
        self.frame.rename(columns=new_names, inplace=True)
        self.frame = self.frame.transpose()
        # Group low abundant into 'others' #
        low_abundance = self.frame.sum() < 30000
        other_count = self.frame.loc[:, low_abundance].sum(axis=1)
        self.frame = self.frame.loc[:, ~low_abundance]
        self.frame['Others'] = other_count
        # Normalize #
        self.frame = self.frame.apply(lambda x: 100*x/x.sum(), axis=1)
        # Sort the table by sum #
        sums = self.frame.sum()
        sums.sort(ascending=False)
        self.frame = self.frame.reindex_axis(sums.keys(), axis=1)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', stacked=True, color=cool_colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_ylabel('Relative abundances in percent')
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.40), fancybox=True, shadow=True, ncol=5, prop={'size':10})
        # Font size #
        axes.tick_params(axis='x', which='major', labelsize=11)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)


################################################################################
# Get the projects #
#proj = illumitag.projects['evaluation']
#proj.load()
proj.graphs += [BarcodeStack(proj)]
proj.graphs += [LengthDistribution(proj)]
proj.graphs += [FractionTaxaBarStack(proj)]
print "proj.graphs[-1].plot()"
print "proj.graphs[-2].plot()"
print "proj.graphs[-3].plot()"

# Get the cluster #
cluster = illumitag.clustering.favorites.evaluation
