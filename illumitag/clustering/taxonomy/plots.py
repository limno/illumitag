# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph

# Third party modules #
from matplotlib import pyplot
import brewer2mpl

# Constants #
__all__ = ['TaxaBarstack']

################################################################################
class TaxaBarstack(Graph):
    """Distribution of named species by sample"""
    short_name = 'taxa_barstack'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Sorting by fraction #
        if self.parent.samples[0].info.get('Filter_fraction'):
            samples = sorted(self.parent.samples, key = lambda s: (s.info['Filter_fraction'], s.short_name))
            self.frame = self.frame.reindex(index=[s.short_name for s in samples])
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
        self.save_plot(fig, axes, width=24.0, height=14.0, bottom=0.20, top=0.97, left=0.04, right=0.98)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
class TaxaHeatmap(Graph):
    """Distribution of named species by sample"""
    short_name = 'taxa_heatmap'

    def plot(self):
        pass