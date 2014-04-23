# Internal modules #
from illumitag.common.autopaths import AutoPaths

# Third party modules #
from rpy2 import robjects as ro

###############################################################################
class BetaDispersion(object):

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
        # Prepare #
        ro.r("library(vegan)")
        ro.r("data = read.table('%s', header=TRUE, sep='\t', row.names='OTUID')" % (self.table.path))
        ro.r("meta = read.table('%s', header=TRUE, sep='\t', row.names=1)" % (self.meta_data_path))
        # Compute #
        ro.r("data_ordered = data[order(row.names(data)),]")
        ro.r("data_sqrt = sqrt(data_ordered)")
        ro.r("data_wa = wisconsin(data_sqrt)")
        ro.r("data_dist = vegdist(data_wa, method='%s')" % self.dist_method)
        # Meta data #
        ro.r("meta = meta[row.names(data),]")
        ro.r("meta_ordered = meta[order(row.names(meta)),]")
        ro.r("pool = factor(meta[,1])")
        ro.r("chemistry = factor(meta[,3])")
        # Chemistry group #
        ro.r("mod1 = betadisper(data_dist, chemistry)")
        ro.r("test = permutest(mod1, control = permControl(nperm = 1000))")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_permutest, 'w') as handle: handle.write(result)
        ro.r("test = anova(mod1)")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_anova, 'w') as handle: handle.write(result)
        ro.r("pdf(file='%s')" % self.p.beta_dispersion_chemistry_plots)
        ro.r("plot(mod1)")
        ro.r("boxplot(mod1)")
        ro.r("plot(TukeyHSD(mod1))")
        ro.r("dev.off()")
        # Pool group #
        ro.r("mod2 = betadisper(data_dist, pool)")
        ro.r("test = permutest(mod2, control = permControl(nperm = 1000))")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_permutest, 'w') as handle: handle.write(result)
        ro.r("test = anova(mod2)")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_anova, 'w') as handle: handle.write(result)
        ro.r("pdf(file='%s')" % self.p.beta_dispersion_pool_plots)
        ro.r("plot(mod2)")
        ro.r("boxplot(mod2)")
        ro.r("plot(TukeyHSD(mod2))")
        ro.r("dev.off()")