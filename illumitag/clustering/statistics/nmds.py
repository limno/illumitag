# Internal modules #
from illumitag.common.tmpstuff import TmpFile
from illumitag.common.autopaths import AutoPaths

# Third party modules #
import sh

###############################################################################
class NMDS(object):

    all_paths = """
    /nmds/
    /permanova/
    /betadis/
    """

    def __init__(self, parent):
        # Save parent #
        self.stat, self.parent = parent, parent
        self.otu = parent.otu
        # Paths #
        self.p = AutoPaths(self.parent.p.betadis_dir, self.all_paths)

    def run(self):
        # Script #
        script = []
        # Load libs #
        script += ["library(vegan)"]
        script += ["library(MASS)"]
        script += ["library(ggplot2)"]
        script += ["library(compare)"]
        # Load data #
        script += ["data = read.table('%s', header=TRUE, sep='\\t', row.names='OTUID')" % (self.table.path)]
        script += ["meta = read.table('%s', header=TRUE, sep='\\t', row.names=1)" % (self.meta_data_path)]
        # Compute nmds #
        script += ["ord = metaMDS(data, distance='%s', trymax=200)" % self.dist_method]
        script += ["nmds = scores(ord)"]
        # Make a dataframe #
        script += ["df = merge(meta, nmds, by.x='row.names', by.y='row.names')"]
        # Make factors #
        script += ["df$barcode = factor(df$barcode)"]
        script += ["df$pool = factor(df$pool)"]
        script += ["df$chemistry = factor(df$chemistry)"]
        # Make plots #
        script += ["p = ggplot(df, aes(NMDS1, NMDS2)) + xlab('Dimension 1') + ylab('Dimension 2')"]
        script += ["pdf(file='%s')" % self.p.nmds_by_barcode]
        script += ["p + geom_point(aes(colour=barcode))"]
        script += ["dev.off()"]
        script += ["pdf(file='%s')" % self.p.nmds_by_pool]
        script += ["p + geom_point(aes(colour=pool))"]
        script += ["dev.off()"]
        script += ["pdf(file='%s')" % self.p.nmds_by_chemistry]
        script += ["p + geom_point(aes(colour=chemistry))"]
        script += ["dev.off()"]
        # Run it #
        sh.R('--no-save', '-f', TmpFile.from_string('\n'.join(script)), _out=self.p.nmds_out)