# Internal modules #
from illumitag.graphs import Graph

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['UniqueStarts']

################################################################################
class UniqueStarts(Graph):
    """When guessing barcodes: how many unique starting kmers ?"""
    short_name = 'unique_starts'
    kmer = 7
    upto = 5000

    def plot(self):
        # Data #
        gen = iter(self.parent)
        self.seqs = [gen.next() for x in xrange(self.upto)]
        self.func = lambda x: len(set([r.seq.tostring()[0:self.kmer] for r in self.seqs[0:x]]))
        self.data = [self.func(x) for x in xrange(self.upto)]
        self.frame = pandas.Series(self.data)
        # Plot #
        fig = pyplot.figure()
        pyplot.plot(self.data, 'k-')
        title = 'Rarefaction curve of unique %i mers.' % self.kmer
        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Number of sequences considered')
        axes.set_ylabel('Unique %i mers' % self.kmer)
        axes.yaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        self.frame.to_csv(self.csv_path)
        pyplot.close(fig)
