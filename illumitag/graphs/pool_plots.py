# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from illumitag.common import flatten
from illumitag.graphs import Graph
from illumitag.helper.chimeras import UchimeRef, UchimeDenovo

# Third party modules #
import matplotlib, pandas
from matplotlib import pyplot
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
__import__('statsmodels.api')
from statsmodels.formula.api import ols
import statsmodels.graphics as smgraphics

# Constants #
__all__ = ['BarcodeHist', 'ReadsThatPassHist', 'SalvageHist', 'MissmatchReg', 'PrimerCounts', 'ReadsWithN', 'QualityFilter', 'LenFilter', 'UchimeRefBar', 'UchimeDenovoBar']

################################################################################
class BarcodeHist(Graph):
    """Distribution of barcodes within a pool"""
    short_name = 'barcode_hist'

    def plot(self):
        # Data #
        self.frame = pandas.Series(self.parent.good_barcodes.breakdown.values(), index=self.parent.samples.bar_names)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', color='gray')
        fig = pyplot.gcf()
        axes.set_title('Distribution of reads with matching barcodes for pool %i by barcode' % self.parent.num)
        axes.set_ylabel('Number of paired reads')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        self.frame.to_csv(self.csv_path)

################################################################################
class ReadsThatPassHist(Graph):
    """Lorem"""
    short_name = 'reads_that_pass'

    def plot(self):
        # Data #
        columns = ["After barcode check", "After processing"]
        rows = self.parent.samples.bar_names
        data = zip(self.parent.good_barcodes.breakdown.values(), self.parent.quality_reads.untrimmed.good_barcodes_breakdown.values())
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar')
        fig = pyplot.gcf()
        axes.set_title('Distribution of barcodes before and after processing for pool %i' % self.parent.num)
        axes.set_ylabel('Number of paired reads')
        axes.yaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, sep=('y'), bottom=0.12)
        self.frame.to_csv(self.csv_path)

################################################################################
class SalvageHist(Graph):
    """Breakdown of outcomes of barcodes within a pool"""
    short_name = 'salvage_hist'

    def plot(self):
        # Data #
        columns = [g.doc for g in self.parent]
        rows = self.parent.samples.bar_sided_names
        data = [[g.counter[sided_name] for g in self.parent.outcomes] for sided_name in rows]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='bar', stacked=True, color=['g','k','y','orange','r'])
        fig = pyplot.gcf()
        axes.set_title('Distribution of barcodes by category for pool %i' % self.parent.num)
        axes.set_ylabel('Number of single reads')
        axes.yaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, sep=('y'), bottom=0.12)
        self.frame.to_csv(self.csv_path)

################################################################################
class MissmatchReg(Graph):
    """Regression of mismatches per barcode proportions"""
    short_name = 'missmatch_reg'

    def plot(self):
        # Data #
        matched = []
        missmatched = []
        names = []
        for first_bar, second_bar in self.parent.samples.all_bar_pairs:
            if frozenset((first_bar, second_bar)) not in self.parent.bad_barcodes.set_counter: continue
            names.append(first_bar[7:] + '-' + second_bar[7:])
            matched.append(self.parent.good_barcodes.counter[first_bar] / len(self.parent.good_barcodes))
            amount_miss = self.parent.bad_barcodes.set_counter[frozenset((first_bar, second_bar))]
            total_miss = self.parent.bad_barcodes.counter[second_bar]
            missmatched.append(amount_miss / total_miss)
        # Regression #
        self.reg = ols("data ~ x", data=dict(data=missmatched, x=matched)).fit()
        # Plot #
        fig = smgraphics.regressionplots.plot_fit(self.reg, 1)
        axes = pyplot.gca()
        smgraphics.regressionplots.abline_plot(model_results=self.reg, ax=axes, color='red')
        axes.set_title('Barcode abundance against mismatched barcode abundance for every barcode pair')
        axes.set_xlabel('Fraction of matched barcode against all matched')
        axes.set_ylabel('Fraction of mismatched barcode against all mismatched')
        axes.legend(loc=1)
        # Add regression result #
        matplotlib.rc('text', usetex=True)
        text = "\\centerline{\\textbf{Least squares fit}}\n"
        text += "{\\raggedright \\textit{Residues}: \\texttt{%i}}\n"
        text += "{\\raggedright \\textit{R\\textsuperscript{2}:} \\texttt{%.2f}}\n"
        text += "{\\raggedright \\textit{Slope}: \\texttt{%.2f}}"
        text = text % (self.reg.df_resid, self.reg.rsquared, self.reg.params[1])
        anchor = AnchoredText(text, prop=dict(size=14), frameon=True, loc=2)
        anchor.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        axes.add_artist(anchor)
        # Find outliers #
        outliers = (idx for idx,t in enumerate(self.reg.outlier_test().icol(2)) if t < 0.5)
        outliers = [(matched[i], missmatched[i], names[i]) for i in outliers]
        for x,y,name in outliers: axes.annotate(name, xy = (x,y), size=9, xytext = (-5, 5),
                                                textcoords = 'offset points', ha = 'center', va = 'bottom',
                                                bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3))
        # Save it #
        self.save_plot(fig, axes)
        matplotlib.rc('text', usetex=False)

################################################################################
class AssemblyCounts(Graph):
    """How many reads assembled and how many did not"""
    short_name = 'assembly_counts'

    def plot(self):
        # Data #
        rows = [bg.doc for bg in self.parent]
        columns = ['Assembled', 'Unassembled', 'Low quality']
        percentage = lambda x,y: (len(x)/len(y))*100 if len(y) != 0 else 0
        data = [(percentage(bg.assembled, bg),
                 percentage(bg.unassembled, bg),
                 (bg.assembled.stats['lowqual']/len(bg))*100) for bg in self.parent]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh', stacked=True)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Proportion of reads pairs that assemble for pool %i' % self.parent.num)
        axes.set_xlabel('Percentage of total reads for a given category of barcode outcome')
        axes.xaxis.grid(True)
        axes.set_xlim(0,100)
        # Save it #
        self.save_plot(fig, axes, left=0.15, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class PrimerCounts(Graph):
    """Are the primers there on the assembled reads"""
    short_name = 'primer_counts'

    def plot(self):
        # Data #
        rows = flatten([(bg.doc + '\n and they assembled', bg.doc + '\n and they are unassembled') for bg in self.parent])
        columns = [pg.__doc__ for pg in self.parent.good_barcodes.assembled.children]
        data = flatten([([len(pg) for pg in bg.assembled],[len(pg) for pg in bg.unassembled]) for bg in self.parent])
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        colors = ['g','r','y','orange','k']
        axes = self.frame.plot(kind='barh', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Primer presence check results for pool %i (0 base pairs missmatches are allowed)' % self.parent.num)
        axes.set_xlabel('Number of paired reads')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, left=0.15, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class ReadsWithN(Graph):
    """How many reads where lost because of N bases"""
    short_name = 'reads_with_n'

    def plot(self):
        # Data #
        rows = flatten([(bg.doc + '\n and they assembled', bg.doc + '\n and they are unassembled') for bg in self.parent.outcomes])
        columns = [pg.__doc__ for pg in self.parent.good_barcodes.assembled.children]
        percentage = lambda x,y: 100-(len(x)/len(y))*100 if len(y) != 0 else 0
        data_ass = [[percentage(pg.n_filtered, pg.orig_reads) for pg in bg.assembled] for bg in self.parent.outcomes]
        data_unass = [[percentage(pg.n_filtered, pg.orig_reads) for pg in bg.unassembled] for bg in self.parent.outcomes]
        data = flatten(zip(data_ass,data_unass))
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        colors = ['g','r','y','orange','k']
        axes = self.frame.plot(kind='barh', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Fraction of reads discarded because of undetermined base pairs for pool %i' % self.parent.num)
        axes.set_xlabel('Percentage of paired reads with "N" (stacked)')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, left=0.15, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class QualityFilter(Graph):
    """How many reads where lost because of low quality"""
    short_name = 'quality_filter'

    def plot(self):
        # Data #
        rows = [bg.doc + '\n and they assembled' for bg in self.parent]
        columns = [pg.__doc__ for pg in self.parent.good_barcodes.assembled.children]
        percentage = lambda x,y: 100-(len(x)/len(y))*100 if len(y) != 0 else 0
        data = [[percentage(pg.qual_filtered, pg.n_filtered) for pg in bg.assembled] for bg in self.parent]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        colors = ['g','r','y','orange','k']
        axes = self.frame.plot(kind='barh', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title("Fraction of assembled reads discarded because of low quality for pool %i (after discarding N's)" % self.parent.num)
        axes.set_xlabel('Percentage of assembled reads with bad PHRED scores (stacked)')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, left=0.15, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class LenFilter(Graph):
    """How many reads where lost because of too much overlap"""
    short_name = 'len_filter'

    def plot(self):
        # Data #
        rows = [bg.doc + '\n and they assembled' for bg in self.parent]
        columns = [pg.__doc__ for pg in self.parent.good_barcodes.assembled.children]
        percentage = lambda x,y: 100-(len(x)/len(y))*100 if len(y) != 0 else 0
        data = [[percentage(pg.len_filtered, pg.qual_filtered) for pg in bg.assembled] for bg in self.parent]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        colors = ['g','r','y','orange','k']
        axes = self.frame.plot(kind='barh', stacked=True, color=colors)
        fig = pyplot.gcf()
        # Other #
        axes.set_title("Fraction of assembled reads discarded because of too much overlap for pool %i (after quality filter)" % self.parent.num)
        axes.set_xlabel('Percentage of assembled reads that were too short (stacked)')
        axes.xaxis.grid(True)
        # Save it #
        self.save_plot(fig, axes, left=0.15, sep=('x'))
        self.frame.to_csv(self.csv_path)

################################################################################
class UchimeRefBar(Graph):
    """Lorem"""
    short_name = "uchime_ref"

    def plot(self):
        # Data #
        rows = [bg.doc + '\n and they assembled\n and the primer was found' for bg in self.parent.outcomes]
        data = [bg.assembled.good_primers.uchime_ref.percent for bg in self.parent.outcomes]
        self.frame = pandas.Series(data, index=rows)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh')
        fig = pyplot.gcf()
        # Other #
        axes.set_title(UchimeRef.title % self.parent.num)
        fig.suptitle("Downsampled to %i" % UchimeRef.downto)
        axes.set_xlabel('Percentage of sequences identified as chimeras after quality filtering')
        axes.yaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, left=0.15)
        self.frame.to_csv(self.csv_path)

################################################################################
class UchimeDenovoBar(Graph):
    """Lorem"""
    short_name = "uchime_denovo"

    def plot(self):
        # Data #
        rows = [bg.doc + '\n and they assembled\n and the primer was found' for bg in self.parent.outcomes]
        data = [bg.assembled.good_primers.uchime_denovo.percent for bg in self.parent.outcomes]
        self.frame = pandas.Series(data, index=rows)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh')
        fig = pyplot.gcf()
        # Other #
        axes.set_title(UchimeDenovo.title % self.parent.num)
        fig.suptitle("Downsampled to %i" % UchimeDenovo.downto)
        axes.set_xlabel('Percentage of sequences identified as chimeras after quality filtering')
        axes.yaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, left=0.15)
        self.frame.to_csv(self.csv_path)
