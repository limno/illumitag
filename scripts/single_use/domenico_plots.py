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
__all__ = ['MainRiver', 'Tributaries']

################################################################################
class MainRiver(Graph):
    """Without the tributaries"""
    short_name = 'main_river'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Drop tributaries #
        key = lambda s: (s.info['Filter_fraction'], s.short_name)
        samples = sorted([s for s in self.parent.samples if s.info['Tributary']=='2'], key=key)
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
class Tributaries(Graph):
    """Only the tributaries"""
    short_name = 'tributaries'

    def plot(self):
        # Data #
        self.frame = self.parent.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
        # Drop the river #
        self.frame = self.frame.loc[[s.short_name for s in self.parent.samples if s.info['Tributary']=='1']]
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
samples = [s for pool in illumitag.runs[4][0:3] for s in pool.samples if s.used]
samples += illumitag.runs[4][3][0:13]
cluster = illumitag.clustering.Cluster(samples, 'domenico')
phyla = cluster.otu_uparse.taxonomy_silva.comp_phyla
river = MainRiver(phyla)
tributaries = Tributaries(phyla)
phyla.graphs += [river]
phyla.graphs += [tributaries]
