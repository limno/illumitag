# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #
import illumitag
from illumitag.graphs import Graph
from illumitag.helper.freshwater import species_names, clade_names

# Third party modules #
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import brewer2mpl, numpy, matplotlib, pandas

# Constants #
__all__ = ['MainRiver', 'Tributaries', 'TaxaHeatmap', 'TaxaHeatmapMainRiver', 'TaxaHeatmapTributaries']

################################################################################
class TaxaBarstack(Graph):
    """Distribution of named species by sample"""
    short_name = 'taxa_barstack'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Take only the ones we use #
        self.samples = [s for s in self.parent.samples if s.short_name in self.parent.taxa_table.index]
        # Drop tributaries or main river #
        self.samples = [s for s in self.samples if s.info['Tributary']==self.tributary]
        # Sort #
        key = lambda s: (s.info['Filter_fraction'], s.short_name)
        self.samples = sorted(self.samples, key=key)
        self.frame = self.frame.reindex(index=[s.short_name for s in self.samples])
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
        axes.set_title('Species relative abundances per sample (blasting against "%s" database)' % self.parent.taxonomy.database)
        axes.set_ylabel('Relative abundances in percent')
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=5)
        # Save it #
        self.save_plot(fig, axes, width=24.0, height=14.0, bottom=0.20, top=0.97, left=0.04, right=0.98)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class MainRiver(TaxaBarstack):
    """Without the tributaries"""
    short_name = 'main_river'
    tributary = '2'

################################################################################
class Tributaries(TaxaBarstack):
    """Only the tributaries"""
    short_name = 'tributaries'
    tributary = '1'

################################################################################
class TaxaHeatmap(Graph):
    """Abundance of most abundant species for every sample"""
    short_name = 'taxa_heatmap'
    targets = ['LD12', 'acI-B1', 'acI-A7', 'acI-C2', 'Pnec', 'Luna1-A2', 'Algor', 'Iluma-A1', 'acI-A4', 'acSTL-A1', 'Iluma-C1']

    def plot(self):
        # Data #
        self.orig_frame = self.parent.taxa_table
        self.norm_frame = self.orig_frame.apply(lambda x: x/x.sum(), axis=1) # Try different normalization ?
        # Take only the right tributary kind #
        self.samples = [s for s in self.parent.samples if s.info['Tributary'] == self.tributary]
        # Take only the ones we use #
        self.samples = [s for s in self.samples if s.short_name in self.orig_frame.index]
        # Sorting by fraction #
        self.samples = sorted(self.samples, key = lambda s: (s.info['Filter_fraction'], s.short_name))
        # Filter the frame #
        sample_names = [s.short_name for s in self.samples]
        self.sort_frame = self.norm_frame.reindex(index=sample_names)
        # Take only our targets and transpose #
        self.targ_frame = self.sort_frame[self.targets]
        self.targ_frame = self.targ_frame.transpose()
        # Plot #
        fig = pyplot.figure()
        fig.set_figwidth(24.0)
        fig.set_figheight(14.0)
        gs = gridspec.GridSpec(2, 2, height_ratios=[1,4], width_ratios=[10,1])
        axes = fig.add_subplot(gs[2])
        heatmap = axes.pcolor(self.targ_frame, cmap=pyplot.cm.Blues, alpha=0.8, edgecolors='#ACABFE')
        # Other #
        axes.grid(False)
        axes.invert_yaxis()
        axes.set_frame_on(False)
        axes.set_yticks(numpy.arange(0.5, len(self.targ_frame.index), 1), minor=False)
        axes.set_xticks(numpy.arange(0.5, len(self.targ_frame.columns), 1), minor=False)
        axes.set_yticklabels(self.targ_frame.index, minor=False)
        axes.set_xticklabels(self.targ_frame.columns, minor=False, rotation='vertical')
        axes.xaxis.tick_bottom()
        axes.yaxis.tick_left()
        axes.get_yaxis().set_tick_params(which='both', direction='out')
        axes.set_title('Species relative abundances per sample (top 12 species-level only)')
        # Extra barstack #
        categories = ['Freshwater taxa', 'Freshwater clade or lineage', 'Everything else']
        colors = [(1.0, 0.49, 0.0), (1.0, 1.0, 0.2), (0.4, 0.76, 0.64)]
        self.breakdown = defaultdict(lambda: defaultdict(int))
        for sample_name, column in self.orig_frame.iterrows():
            for taxa_name, count in column.iteritems():
                if taxa_name in species_names: category = categories[0]
                elif taxa_name in clade_names: category = categories[1]
                else:                          category = categories[2]
                self.breakdown[category][sample_name] += count
        self.breakdown = pandas.DataFrame(self.breakdown)
        self.breakdown = self.breakdown.fillna(0)
        self.breakdown = self.breakdown.apply(lambda x: x/x.sum(), axis=1)
        self.breakdown = self.breakdown.reindex(index=sample_names)
        axes = fig.add_subplot(gs[0], sharex=axes)
        self.breakdown.plot(kind='bar', stacked=True, ax=axes, color=colors)
        percentage = lambda x, pos: '%1.0f%%' % (x*100.0)
        axes.set_ylim([0,1])
        axes.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.patch.set_alpha(0.0)
        axes.spines['right'].set_visible(False)
        axes.spines['top'].set_visible(False)
        axes.spines['bottom'].set_visible(False)
        axes.tick_params(axis='both', direction='out')
        axes.xaxis.tick_bottom()
        axes.yaxis.tick_left()
        axes.set_title('Classification level breakdown per sample (blasting against "%s" database)' % self.parent.taxonomy.database)
        # Extra scale #
        axes = fig.add_subplot(gs[3])
        fig.colorbar(heatmap, pad=0.1, fraction=1.0, shrink=1.0, format=matplotlib.ticker.FuncFormatter(percentage))
        axes.set_axis_off()
        # Save it #
        fig.savefig(self.path)
        self.targ_frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class TaxaHeatmapMainRiver(TaxaHeatmap):
    tributary = '2'
    short_name = 'taxa_heatmap_main_river'

################################################################################
class TaxaHeatmapTributaries(TaxaHeatmap):
    tributary = '1'
    short_name = 'taxa_heatmap_tributaries'

################################################################################
# Make a cluster of samples #
cluster = illumitag.clustering.favorites.danube
phyla = cluster.otu_uparse.taxonomy_silva.comp_phyla
river = MainRiver(phyla)
tributaries = Tributaries(phyla)
phyla.graphs += [river]
phyla.graphs += [tributaries]

tips = cluster.otu_uparse.taxonomy_fw.comp_tips
heatmap_river = TaxaHeatmapMainRiver(tips)
heatmap_tributaries = TaxaHeatmapTributaries(tips)
tips.graphs += [heatmap_river]
tips.graphs += [heatmap_tributaries]

#cluster.otu_uparse.taxonomy_silva.comp_phyla.make_plots()
#cluster.otu_uparse.taxonomy_fw.comp_tips.make_plots()
