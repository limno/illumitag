# Built-in modules #

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import flatten

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['BarcodeStack', 'AssemblyCounts']

################################################################################
class BarcodeStack(Graph):
    """General distribution of barcodes for all pools"""
    short_name = 'barcode_stack'

    def plot(self):
        # Data #
        rows = ['Pool "%s" (%s)' % (p.long_name, p) for p in self.parent.pools]
        columns = [o.doc for o in self.parent.first.outcomes]
        data = [[o.count for o in p.outcomes] for p in self.parent.pools]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh', stacked=True, color=['g','k','y','orange','r'])
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Demultiplexing result (%i barcodes are used)' % len(self.parent.first.samples))
        axes.set_xlabel('Number of paired reads')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class AssemblyCounts(Graph):
    """Do the reads assemble together or not"""
    short_name = 'assembly_counts'

    def plot(self):
        # Data #
        rows = ['Pool "%s" (%s)' % (p.long_name, p) for p in self.parent.pools]
        columns = flatten([(o.short_name + "_ass", o.short_name + "_unass") for o in self.parent.first.outcomes])
        data = [flatten([(len(o.assembled), len(o.unassembled)) for o in pool.outcomes]) for pool in self.parent.pools]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        colors = ['g','g','gray','gray','y','y','orange','orange','r','r']
        axes = self.frame.plot(kind='barh', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Add pattern #
        unass_patches = [p for i,p in enumerate(axes.patches) if (i/len(self.parent))%2 == 1]
        for p in unass_patches: p.set_hatch('//')
        # Other #
        axes.set_title('Assembling reads using PANDAseq')
        axes.set_xlabel('Number of paired reads')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
        self.frame.to_csv(self.csv_path)
