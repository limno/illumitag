# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph, cool_colors

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['TaxaBarstack']

################################################################################
class TaxaBarstack(Graph):
    """Distribution of named species by sample"""
    short_name = 'taxa_barstack'
    width = 24.0
    height = 14.0
    bottom = 0.20
    top = 0.97
    left = 0.04
    right = 0.98
    legend_anchor = -0.10

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', stacked=True, color=cool_colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Species relative abundances per sample (blasting against "%s" database)' % self.parent.taxonomy.database)
        axes.set_ylabel('Relative abundances in percent')
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, self.legend_anchor), fancybox=True, shadow=True, ncol=5)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)