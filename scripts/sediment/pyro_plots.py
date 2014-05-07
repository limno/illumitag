# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
import illumitag
from illumitag.graphs import Graph
from sediment_plots import sediment_new_names

# Third party modules #
from matplotlib import pyplot
import brewer2mpl

# Constants #
__all__ = ['PyroSorted']

################################################################################
class PyroSorted(Graph):
    """Sorted pyro samples"""
    short_name = 'sorted_pyro'
    width = 22.0
    height = 14
    bottom = 0.3
    top = 0.95
    left = 0.05
    right = 0.95
    formats = ('pdf', 'eps')

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Rename #
        self.frame = self.frame.reindex(index=sediment_new_names.keys())
        self.frame = self.frame.rename(index=sediment_new_names)
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
        axes.set_ylabel('Relative abundances in percent')
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        axes.set_ylim([0,100])
        # Put a legend below current axis
        axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), fancybox=True, shadow=True, ncol=5)
        # Save it #
        self.save_plot(fig, axes)
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
# Make a cluster of samples #
cluster = illumitag.clustering.favorites.pyro_comparison

# Get phyla #
phyla = cluster.otu_uparse.taxonomy_silva.comp_phyla
pyro_sorted = PyroSorted(phyla)
phyla.graphs += [pyro_sorted]
print "phyla.graphs[-1].plot()"