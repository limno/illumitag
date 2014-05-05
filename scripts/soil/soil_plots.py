# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
import illumitag
from illumitag.graphs import Graph

# Third party modules #
import pandas
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
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)

################################################################################
# Get the projects #
#proj = illumitag.projects['evaluation']
#proj.load()
proj.graphs += [BarcodeStack(proj)]
#print "proj.graphs[-1].plot()"

# Get the cluster #
cluster = illumitag.clustering.favorites.evaluation
