#!/usr/bin/env python2

"""
A script to make some graphs for the soda project.
"""

# Modules #
from __future__ import division
from matplotlib import pyplot
import illumitag, math, brewer2mpl
from pandas.rpy.common import convert_to_r_dataframe
from rpy2 import robjects as ro

###############################################################################
temporal = [s for s in illumitag.runs[3][7].samples if s.used]
spatial = [s for s in illumitag.runs[4][3].samples if s.used and s.group_name != 'JDS_2007']
samples = temporal + spatial
cluster = illumitag.clustering.Cluster(samples, 'soda')

colors = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
colors.reverse()
colors += brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
colors += brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors

###############################################################################
def temp_cond_to_sal(cond, temp):
    """http://www.fivecreeks.org/monitor/sal.shtml"""
    a0 = 0.008
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 = 2.7081
    b0 = 0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 = 0.0636
    b5 = -0.0144
    c0 = 0.6766097
    c1 = 0.0200564
    c2 = 0.0001104259
    c3 = -0.00000069698
    c4 = 0.0000000010031
    if temp < 0.0 or temp > 30.0: return None #raise Exception("Out of range")
    if cond < 0.0: raise Exception("Out of range")
    r = cond/42914;
    r /= (c0+temp*(c1+temp*(c2+temp*(c3+temp*c4))))
    r2 = math.sqrt(r);
    ds = b0+r2*(b1+r2*(b2+r2*(b3+r2*(b4+r2*b5))));
    ds *= ((temp-15.0)/(1.0+0.0162*(temp-15.0)));
    sal = a0+r2*(a1+r2*(a2+r2*(a3+r2*(a4+r2*a5))))+ds;
    if sal < 2.0: return None #raise Exception("Under scale")
    if sal > 42.0: raise Exception("Over scale")
    return sal

###############################################################################
def temp_graph():
    temps = [s.info['temperature'][0] for s in temporal]
    fig, axes = pyplot.subplots()
    fig.set_figwidth(8)
    fig.set_figheight(1)
    fig.subplots_adjust(bottom=0.3, top=0.98, left=0.02, right=0.98)
    axes.plot(temps, [1]*len(temps), 'ko')
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.xaxis.tick_bottom()
    axes.yaxis.set_visible(False)
    #axes.set_xlabel('Temperature [Celsius]')
    fig.savefig("temp.pdf")

def ph_graph():
    ph = [s.info['pH'][0] for s in samples if 'pH' in s.info and s.info['pH']]
    fig, axes = pyplot.subplots()
    fig.set_figwidth(8)
    fig.set_figheight(1)
    fig.subplots_adjust(bottom=0.3, top=0.98, left=0.02, right=0.98)
    axes.plot(ph, [1]*len(ph), 'ko')
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.xaxis.tick_bottom()
    axes.yaxis.set_visible(False)
    #axes.set_xlabel('pH')
    fig.savefig("ph.pdf")

def salt_graph():
    for s in samples:
        if not s.info.get('conductance') or not s.info.get('temperature'): continue
        salinity = temp_cond_to_sal(s.info['conductance'][0], s.info['temperature'][0])
        print 'cond: %s, temp: %s, sal: %s' % (s.info['conductance'][0], s.info['temperature'][0], salinity)
        if salinity: s.info['salinity'] = salinity
    sals = [s.info['salinity'] for s in samples if 'salinity' in s.info]
    fig, axes = pyplot.subplots()
    fig.set_figwidth(8)
    fig.set_figheight(1)
    fig.subplots_adjust(bottom=0.3, top=0.98, left=0.02, right=0.98)
    axes.plot(sals, [1]*len(sals), 'ko')
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.xaxis.tick_bottom()
    axes.yaxis.set_visible(False)
    #axes.set_xlabel('Salinity')
    fig.savefig("sals.pdf")

def ph_temp_graph():
    with_ph = [s for s in samples if 'pH' in s.info and s.info['pH']]
    x = [s.info['pH'][0] for s in with_ph]
    y = [s.info['temperature'][0] for s in with_ph]
    fig, axes = pyplot.subplots()
    fig.set_figwidth(8)
    fig.set_figheight(8)
    fig.subplots_adjust(bottom=0.3, top=0.98, left=0.06, right=0.98)
    axes.plot(x, y, 'ko')
    axes.set_xlabel('pH')
    axes.set_ylabel('Temp')
    # Add annotations #
    for i in range(len(with_ph)):
        pyplot.annotate(with_ph[i].short_name, size=9, xy = (x[i], y[i]), xytext = (10, 0),
                        textcoords = 'offset points', ha = 'left', va = 'center',
                        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
    # Save #
    fig.savefig("ph_against_temp.pdf")

###############################################################################
def temporal_taxa_barstack():
    # Data #
    frame = cluster.otu_uparse.taxonomy_silva.comp_phyla.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
    # Sorting by fraction #
    samples = sorted(temporal, key = lambda s: (s.info['group'], s.info['num']))
    frame = frame.reindex(index=[s.short_name for s in samples])
    # Plot #
    fig = pyplot.figure()
    axes = frame.plot(kind='bar', stacked=True, color=colors)
    fig = pyplot.gcf()
    # Other #
    axes.set_title('Species relative abundances per sample')
    axes.set_ylabel('Relative abundances in percent')
    axes.xaxis.grid(False)
    axes.yaxis.grid(False)
    axes.set_ylim([0,100])
    # Put a legend below current axis
    axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=5)
    # Save it #
    fig.set_figwidth(16)
    fig.set_figheight(8)
    fig.subplots_adjust(bottom=0.30, top=0.97, left=0.04, right=0.98)
    fig.savefig("temporal_taxa_barstack.pdf")
    pyplot.close(fig)

###############################################################################
def spatial_taxa_barstack():
    # Data #
    frame = cluster.otu_uparse.taxonomy_silva.comp_phyla.taxa_table.apply(lambda x: 100*x/x.sum(), axis=1)
    # Sorting by fraction #
    march = [s for s in spatial if not 'rerun' in s.info]
    march = sorted(march, key = lambda s: (s.info['pH']))
    frame = frame.reindex(index=[s.short_name for s in march])
    # Plot #
    fig = pyplot.figure()
    axes = frame.plot(kind='bar', stacked=True, color=colors)
    fig = pyplot.gcf()
    # Other #
    axes.set_title('Species relative abundances per sample')
    axes.set_ylabel('Relative abundances in percent')
    axes.xaxis.grid(False)
    axes.yaxis.grid(False)
    axes.set_ylim([0,100])
    # Put a legend below current axis
    axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=5)
    # Save it #
    fig.set_figwidth(16)
    fig.set_figheight(8)
    fig.subplots_adjust(bottom=0.30, top=0.97, left=0.04, right=0.98)
    fig.savefig("spatial_taxa_barstack.pdf")
    pyplot.close(fig)

###############################################################################
def richness():
    # Taxa table #
    otu_table = cluster.otu_uparse.taxonomy_silva.otu_table
    # Compute richness #
    #ro.r.library("vegan")
    #R_otu_table = convert_to_r_dataframe(otu_table)
    #R_result = ro.r.qvalue(R_frame.rx2('pvalues'))
    #qvalues = list(R_result.rx2('qvalues'))
