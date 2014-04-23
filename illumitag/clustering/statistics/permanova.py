# Built-in modules #
import itertools

# Internal modules #
from illumitag.common.autopaths import AutoPaths

# Third party modules #
from rpy2 import robjects as ro

###############################################################################
class PERMANOVA(object):

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    """

    def __init__(self, parent):
        # Save parent #
        self.stat, self.parent = parent, parent
        self.tax = parent.tax
        # Paths #
        self.p = AutoPaths(self.parent.p.betadis_dir, self.all_paths)

    def run(self):
        # Basic PERMANOVA #
        ro.r("library(vegan)")
        ro.r("data = read.table('%s', header=TRUE, sep='\t', row.names='OTUID')" % (self.table.path))
        ro.r("meta = read.table('%s', header=TRUE, sep='\t', row.names=1)" % (self.meta_data_path))
        ro.r("data_ordered = data[order(row.names(data)),]")
        # Meta data #
        ro.r("meta = meta[row.names(data),]")
        ro.r("meta_ordered = meta[order(row.names(meta)),]")
        ro.r("meta_ordered$pool = factor(meta_ordered$pool)")
        ro.r("meta_ordered$barcode = factor(meta_ordered$barcode)")
        ro.r("meta_ordered$chemistry = factor(meta_ordered$chemistry)")
        # Run test #
        ro.r("permanova = adonis(formula = data ~ pool * barcode, data=meta_ordered, permutations=1000, method='%s')" % self.dist_method)
        # Get output #
        result = '\n'.join(ro.r("capture.output(print(permanova))")).encode('utf-8')
        with open(self.p.permanova_all, 'w') as handle: handle.write(result)
        # All pool pairs #
        with open(self.p.permanova_pairs, 'w') as handle:
            for pair in itertools.combinations(['1', '2', '3', '4', '5'], 2):
                ro.r("subset = row.names(meta_ordered[meta_ordered$pool == 1 | meta_ordered$pool == 3,])")
                ro.r("subdata = data_ordered[subset,]")
                ro.r("submeta = meta_ordered[subset,]")
                ro.r("permanova = adonis(formula = subdata ~ pool * barcode, data=submeta, permutations=1000, method='%s')" % self.dist_method)
                handle.write("\n\n ---- Pool %s against Pool %s ---- \n\n" % pair)
                result = '\n'.join(ro.r("capture.output(print(permanova))")).encode('utf-8')
                handle.write(result)