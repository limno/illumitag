# Built-in modules #
import shutil, json

# Internal modules #
from util import save_plot
from common import flatten

# Third party modules #
import pandas
from matplotlib import pyplot

################################################################################
def assembly_distrib(group):
    """Distribution of assembly lengths"""
    # Data #
    values = pandas.Series(group.assembled.stats['lengths'])
    # Plot #
    fig = pyplot.figure()
    axes = values.hist(color='gray', bins=max(values))
    fig = pyplot.gcf()
    title = 'Distribution of sequence overlap between pairs that assemble from PANDAseq output'
    title += ' (pool %i, group "%s")' % (group.pool.num, group.short_name)
    axes.set_title(title)
    axes.set_xlabel('Number of nucleotides in overlap')
    axes.set_ylabel('Number of sequence pairs with this much overlap')
    axes.xaxis.grid(False)
    # Change ticks #
    import matplotlib.ticker as mticker
    myLocator = mticker.MultipleLocator(10)
    axes.xaxis.set_major_locator(myLocator)
    # Save it #
    save_plot(fig, axes, group.p.assembly_distrib_pdf, sep=('y'))
    values.to_csv(group.p.assembly_distrib_csv)
    # Copy length pictures #
    destination = group.pool.experiment.p.assembly_length_2_dir
    destination += "pool%i" % group.pool.num + "_" + group.short_name + '.pdf'
    shutil.copy(group.p.assembly_distrib_pdf, destination)

################################################################################
def ass_primer_pos(group):
    """Distribution of where we find the primer positions
    in the assembled reads"""
    # Data #
    fwd_pos, rev_pos = group.assembled.primer_positions
    if not fwd_pos and not rev_pos: return
    fwd_data = flatten([[k]*v for k,v in fwd_pos.items()])
    rev_data = flatten([[k]*v for k,v in rev_pos.items()])
    # Plot #
    fig = pyplot.figure()
    bins = range(min(rev_data + [0]), max(fwd_data + [0])+1, 1)
    if fwd_pos: pyplot.hist(fwd_data, bins=bins, histtype='stepfilled', color='b', alpha=0.5, label='Forward')
    if rev_pos: pyplot.hist(rev_data, bins=bins, histtype='stepfilled', color='r', alpha=0.5, label='Reverse')
    title = "Distribution of primer positions within assembled sequences"
    title += ' (pool %i, group "%s")' % (group.pool.num, group.short_name)
    axes = pyplot.gca()
    axes.set_title(title)
    axes.set_xlabel('Relative position at which the primer is found')
    axes.set_ylabel('Number of primers found at this position')
    axes.xaxis.grid(False)
    axes.legend()
    # Save it #
    save_plot(fig, axes, group.p.primers_assembled_pdf, sep=('y'))
    json.dump({'fwd':fwd_pos, 'rev':rev_pos}, open(group.p.primers_assembled_json, 'w'))
    # Copy to aggregate dir #
    destination = group.pool.experiment.p.primer_position_dir
    destination += "pool%i" % group.pool.num + "_" + group.short_name + '_ass.pdf'
    shutil.copy(group.p.primers_assembled_pdf, destination)

################################################################################
def unass_primer_pos(group):
    """Distribution of where we find the primer positions
    in the unassembled reads"""
    # Data #
    fwd_pos, rev_pos = group.unassembled.primer_positions
    if not fwd_pos and not rev_pos: return
    fwd_data = flatten([[k]*v for k,v in fwd_pos.items()])
    rev_data = flatten([[k]*v for k,v in rev_pos.items()])
    # Plot #
    fig = pyplot.figure()
    bins = range(min(rev_data + [0]), max(fwd_data + [0])+1, 1)
    if fwd_pos: pyplot.hist(fwd_data, bins=bins, histtype='stepfilled', color='b', alpha=0.5, label='Forward')
    if rev_pos: pyplot.hist(rev_data, bins=bins, histtype='stepfilled', color='r', alpha=0.5, label='Reverse')
    title = "Distribution of primer positions within unassembled sequences"
    title += ' (pool %i, group "%s")' % (group.pool.num, group.short_name)
    axes = pyplot.gca()
    axes.set_title(title)
    axes.set_xlabel('Relative position at which the primer is found')
    axes.set_ylabel('Number of primers found at this position')
    axes.xaxis.grid(False)
    axes.legend()
    # Save it #
    save_plot(fig, axes, group.p.primers_unassembled_pdf, sep=('y'))
    json.dump({'fwd':fwd_pos, 'rev':rev_pos}, open(group.p.primers_unassembled_json, 'w'))
    # Copy to aggregate dir #
    destination = group.pool.experiment.p.primer_position_dir
    destination += "pool%i" % group.pool.num + "_" + group.short_name + '_unass.pdf'
    shutil.copy(group.p.primers_unassembled_pdf, destination)
