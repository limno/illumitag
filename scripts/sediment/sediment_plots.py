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
from scipy.cluster.hierarchy import linkage, dendrogram

# Constants #
sediment_new_names = OrderedDict([
    ("SS15",          "Soda ILL SS15"),
    ("soda_454_SS15", "Soda 454 SS15"),
    ("SS16",          "Soda ILL SS16"),
    ("soda_454_SS16", "Soda 454 SS16"),
    ("SS17",          "Soda ILL SS17"),
    ("soda_454_SS17", "Soda 454 SS17"),
    ("US15",          "Soda ILL US15"),
    ("soda_454_US15", "Soda 454 US15"),
    ("US16",          "Soda ILL US16"),
    ("soda_454_US16", "Soda 454 US16"),
    ("US17",          "Soda ILL US17"),
    ("soda_454_US17", "Soda 454 US17"),
    ("ZL15",          "Soda ILL ZL15"),
    ("soda_454_ZL15", "Soda 454 ZL15"),
    ("ZL16",          "Soda ILL ZL16"),
    ("soda_454_ZL16", "Soda 454 ZL16"),
    ("ZL17",          "Soda ILL ZL17"),
    ("soda_454_ZL17", "Soda 454 ZL17"),
    ("p1bc01",        "Sediment ILL 01"),
    ("p1bc02",        "Sediment ILL 02"),
    ("p1bc03",        "Sediment ILL 03"),
    ("p1bc04",        "Sediment ILL 04"),
    ("p1bc05",        "Sediment ILL 05"),
    ("p1bc06",        "Sediment ILL 06"),
    ("p1bc07",        "Sediment ILL 07"),
    ("p1bc08",        "Sediment ILL 08"),
    ("p2bc01",        "Sediment ILL 09"),
    ("p2bc02",        "Sediment ILL 10"),
    ("p2bc03",        "Sediment ILL 11"),
    ("p2bc04",        "Sediment ILL 12"),
    ("p2bc05",        "Sediment ILL 13"),
    ("p2bc06",        "Sediment ILL 14"),
    ("p2bc07",        "Sediment ILL 15"),
    ("p2bc08",        "Sediment ILL 16"),
    ("p5bc01",        "Sediment ILL 17"),
    ("p5bc02",        "Sediment ILL 18"),
    ("p5bc03",        "Sediment ILL 19"),
    ("p5bc04",        "Sediment ILL 20"),
    ("p5bc05",        "Sediment ILL 21"),
    ("p5bc06",        "Sediment ILL 22"),
    ("p5bc07",        "Sediment ILL 23"),
    ("p5bc08",        "Sediment ILL 24"),
    ("sediment_454",  "Sediment 454 00")
])

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
        pyplot.bar(counts.keys(), counts.values(), 1.0, color='gray', align='center', label='Reads from sediment sample')
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
class Dendogram(Graph):
    """Compare the pyro and illumina samples in a tree using Unifrac distance"""
    short_name = 'unifrac_dendogram'
    width = 10.0
    height = 3.5
    bottom = 0.4
    top = 0.95
    left = 0.1
    right = 0.95
    formats = ('pdf', 'eps')

    def plot(self):
        # The redundant distance matrix #
        self.frame = self.parent.distances.copy()
        assert (self.frame.index==self.frame.columns).all()
        # Rename #
        self.frame.rename(columns=sediment_new_names, index=sediment_new_names, inplace=True)
        # Hierarchical clustering UPGMA #
        clusters = linkage(self.frame, method='average')
        dendrogram(clusters, labels=self.frame.index)
        axes = pyplot.gca()
        fig = pyplot.gcf()
        # Other #
        axes.set_ylabel('Cophenetic distance')
        # Change font #
        mono = matplotlib.font_manager.FontProperties(family='monospace', size=9)
        for label in axes.get_xticklabels(): label.set_fontproperties(mono)
        # Save it #
        self.save_plot(fig, axes)
        pyplot.close(fig)

################################################################################
class UnifracNMDS(Graph):
    short_name = 'nmds'
    width = 8.0
    height = 8.0
    bottom = 0.10
    top = 0.95
    left = 0.12
    right = 0.95
    formats = ('pdf', 'eps')

    def plot(self):
        # Coord #
        x = self.parent.coords['NMDS1'].values
        y = self.parent.coords['NMDS2'].values
        names = self.parent.coords['NMDS1'].keys()
        # Make scatter #
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.set_xlabel('Dimension 1')
        axes.set_ylabel('Dimension 2')
        # Plot every point #
        for i, name in enumerate(names):
            x = self.parent.coords['NMDS1'][name]
            y = self.parent.coords['NMDS2'][name]
            label = sediment_new_names[name]
            xytext = (10, 0)
            ha = 'left'
            if label.startswith('Sediment ILL'):
                axes.plot(x, y, 'ro')
                continue
            if label.startswith('Sediment 454'):
                axes.plot(x, y, 'bo')
                label = 'Sediment 454'
            if label.startswith('Soda ILL'):
                axes.plot(x, y, 'o', color='darkgreen')
            if label.startswith('Soda 454'):
                axes.plot(x, y, 'o', color='lightgreen')
                xytext = (-10, 0)
                ha = 'right'
            pyplot.annotate(label, size=9, xy = (x, y), xytext = xytext,
                            textcoords = 'offset points', ha = ha, va = 'center',
                            bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)
        pyplot.close(fig)

################################################################################
if __name__ == '__main__':
    # Get the projects #
    #proj = illumitag.projects['evaluation']
    #proj.load()
    #proj.graphs += [BarcodeStack(proj)]
    #proj.graphs += [LengthDistribution(proj)]
    #proj.graphs += [FractionTaxaBarStack(proj)]
    #print "proj.graphs[-1].plot()"
    #print "proj.graphs[-2].plot()"
    #print "proj.graphs[-3].plot()"

    # Get the cluster #
    cluster = illumitag.clustering.favorites.pyro_comparison
    unifrac = cluster.otu_uparse.taxonomy_silva.stats.unifrac
    unifrac.tree = Dendogram(unifrac)
    print "unifrac.tree.plot()"

    # Supplementary NMDS #
    unifrac.nmds.graph = UnifracNMDS(unifrac.nmds, base_dir=unifrac.nmds.base_dir)
    print "unifrac.nmds.run()"
    print "unifrac.nmds.graph.plot()"