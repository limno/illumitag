# Built-in modules #
import json

# Internal modules #
from illumitag.graphs import Graph
from illumitag.common import flatten

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['AssemblyDistrib', 'AssemblyPrimerPos', 'UnassemblyPrimerPos']

################################################################################
class AssemblyDistrib(Graph):
    """Distribution of assembly lengths"""
    short_name = 'assembly_distrib'

    def plot(self):
        # Data #
        self.frame = pandas.Series(self.parent.assembled.stats['lengths'])
        # Plot #
        pyplot.figure()
        axes = self.frame.hist(color='gray', bins=max(self.frame))
        fig = pyplot.gcf()
        title = 'Distribution of sequence overlap between pairs that assemble from PANDAseq output'
        title += ' (pool %i, group "%s")' % (self.parent.pool.num, self.parent.short_name)
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in overlap')
        axes.set_ylabel('Number of sequence pairs with this much overlap')
        axes.xaxis.grid(False)
        # Change ticks #
        import matplotlib.ticker as mticker
        myLocator = mticker.MultipleLocator(10)
        axes.xaxis.set_major_locator(myLocator)
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        self.frame.to_csv(self.csv_path)

################################################################################
class AssemblyPrimerPos(Graph):
    """Distribution of where we find the primer positions
    in the assembled reads"""
    short_name = 'ass_primer_pos'

    def plot(self):
        # Data #
        fwd_pos, rev_pos = self.parent.assembled.primer_positions
        if not fwd_pos and not rev_pos: return
        fwd_data = flatten([[k]*v for k,v in fwd_pos.items()])
        rev_data = flatten([[k]*v for k,v in rev_pos.items()])
        # Plot #
        fig = pyplot.figure()
        bins = range(min(rev_data + [0]), max(fwd_data + [0])+1, 1)
        if fwd_pos: pyplot.hist(fwd_data, bins=bins, histtype='stepfilled', color='b', alpha=0.5, label='Forward')
        if rev_pos: pyplot.hist(rev_data, bins=bins, histtype='stepfilled', color='r', alpha=0.5, label='Reverse')
        title = "Distribution of primer positions within assembled sequences"
        title += ' (pool %i, group "%s")' % (self.parent.pool.num, self.parent.short_name)
        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Relative position at which the primer is found')
        axes.set_ylabel('Number of primers found at this position')
        axes.xaxis.grid(False)
        axes.legend()
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        json.dump({'fwd':fwd_pos, 'rev':rev_pos}, open(self.json_path, 'w'))

################################################################################
class UnassemblyPrimerPos(Graph):
    """Distribution of where we find the primer positions
    in the unassembled reads"""
    short_name = 'unass_primer_pos'

    def plot(self):
        # Data #
        fwd_pos, rev_pos = self.parent.unassembled.primer_positions
        if not fwd_pos and not rev_pos: return
        fwd_data = flatten([[k]*v for k,v in fwd_pos.items()])
        rev_data = flatten([[k]*v for k,v in rev_pos.items()])
        # Plot #
        fig = pyplot.figure()
        bins = range(min(rev_data + [0]), max(fwd_data + [0])+1, 1)
        if fwd_pos: pyplot.hist(fwd_data, bins=bins, histtype='stepfilled', color='b', alpha=0.5, label='Forward')
        if rev_pos: pyplot.hist(rev_data, bins=bins, histtype='stepfilled', color='r', alpha=0.5, label='Reverse')
        title = "Distribution of primer positions within unassembled sequences"
        title += ' (pool %i, group "%s")' % (self.parent.pool.num, self.parent.short_name)
        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Relative position at which the primer is found')
        axes.set_ylabel('Number of primers found at this position')
        axes.xaxis.grid(False)
        axes.legend()
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        json.dump({'fwd':fwd_pos, 'rev':rev_pos}, open(self.json_path, 'w'))
