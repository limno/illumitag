# Built-in modules #
import time, getpass, locale

# Internal modules #
from util import save_plot
from common import flatten

# Third party modules #
import matplotlib, pandas
from matplotlib import pyplot
__import__('statsmodels.api')

################################################################################
def barcode_barstack(experiment):
    """General distribution of barcodes for all pools"""
    # Data #
    rows = ['pool %i' % p.num for p in experiment.pools]
    columns = [g.doc for g in experiment.first.groups]
    data = [[g.count for g in p.groups] for p in experiment.pools]
    experiment.barstack = pandas.DataFrame(data, index=rows, columns=columns)
    # Plot #
    fig = pyplot.figure()
    axes = experiment.barstack.plot(kind='barh', stacked=True, color=['g','k','y','orange','r'])
    fig = pyplot.gcf()
    # Other #
    axes.set_title('Demultiplexing result (%i barcodes are used)' % len(experiment.first.barcodes))
    axes.set_xlabel('Number of paired reads')
    axes.xaxis.grid(True)
    # Common #
    fig.set_figwidth(18.0)
    fig.set_figheight(10.0)
    fig.subplots_adjust(hspace=0.0, bottom=0.07, top=0.93, left=0.04, right=0.98)
    fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
    fig.text(0.01, 0.98, 'user: ' + getpass.getuser(), horizontalalignment='left')
    locale.setlocale(locale.LC_ALL, '')
    seperate = lambda x,pos: locale.format("%d", x, grouping=True)
    axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
    # Save it #
    fig.savefig(experiment.p.barcode_stats_pdf)
    experiment.barstack.to_csv(experiment.p.barcode_stats_csv)

################################################################################
def assembly_counts(experiment):
    """Do the reads assemble together or not"""
    # Data #
    rows = ['Pool %i' % pool.num for pool in experiment.pools]
    columns = flatten([(g.short_name + "_ass", g.short_name + "_unass") for g in experiment.first.groups])
    data = [flatten([(len(g.assembled), len(g.unassembled)) for g in pool.groups]) for pool in experiment.pools]
    experiment.assembly_counts = pandas.DataFrame(data, index=rows, columns=columns)
    # Plot #
    fig = pyplot.figure()
    colors = ['g','g','gray','gray','y','y','orange','orange','r','r']
    axes = experiment.assembly_counts.plot(kind='barh', stacked=True, color=colors)
    fig = pyplot.gcf()
    # Add pattern #
    unass_patches = [p for i,p in enumerate(axes.patches) if (i/len(experiment.pools))%2 == 1]
    for p in unass_patches: p.set_hatch('//')
    # Other #
    axes.set_title('Assembling reads using PANDAseq')
    axes.set_xlabel('Number of paired reads')
    axes.xaxis.grid(True)
    # Save it #
    save_plot(fig, axes, experiment.p.assembly_counts_pdf, sep=('x'))
    experiment.assembly_counts.to_csv(experiment.p.assembly_counts_csv)
