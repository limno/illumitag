# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
import illumitag
from illumitag.graphs import Graph

# Third party modules #
from matplotlib import pyplot
import brewer2mpl

# Constants #
__all__ = ['PyroSorted']

################################################################################
class PyroSorted(Graph):
    """Sorted pyro samples"""
    short_name = 'sorted_pyro'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        index = [
                 'soda_454_SS15',
                 'SS15',
                 'soda_454_SS16',
                 'SS16',
                 'soda_454_SS17',
                 'SS17',
                 'soda_454_US15',
                 'US15',
                 'soda_454_US16',
                 'US16',
                 'soda_454_US17',
                 'US17',
                 'soda_454_ZL15',
                 'ZL15',
                 'soda_454_ZL16',
                 'ZL16',
                 'soda_454_ZL17',
                 'ZL17',
                 'p1bc01',
                 'p1bc02',
                 'p1bc03',
                 'p1bc04',
                 'p1bc05',
                 'p1bc06',
                 'p1bc07',
                 'p1bc08',
                 'p2bc01',
                 'p2bc02',
                 'p2bc03',
                 'p2bc04',
                 'p2bc05',
                 'p2bc06',
                 'p2bc07',
                 'p2bc08',
                 'p5bc01',
                 'p5bc02',
                 'p5bc03',
                 'p5bc04',
                 'p5bc05',
                 'p5bc06',
                 'p5bc07',
                 'p5bc08',
                 'soil_454',
                 ]
        self.frame = self.frame.reindex(index=index)
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
# Make a cluster of samples #
cluster = illumitag.clustering.favorites.pyro_comparison

# Get phyla #
phyla = cluster.otu_uparse.taxonomy_silva.comp_phyla
pyro_sorted = PyroSorted(phyla)
phyla.graphs += [pyro_sorted]
