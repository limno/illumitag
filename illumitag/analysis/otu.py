# Built-in modules #
import os, csv, random, itertools
from itertools import izip
from collections import Counter

# Internal modules #
from illumitaq.common import AutoPaths, JobRunner, TmpFile
from illumitaq.common import flatten
from illumitaq.otu_plots import cluster_distribution, otu_distribution, sample_sums, otu_sums
from illumitaq.otu_plots import taxa_hist

# Third party modules #
import shutil, sh, pandas
from rpy2 import robjects as ro

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class OTUs(JobRunner):
    dist_method = 'horn'

    default_steps = [
        #{'pick_otus':                 {}},
        #{'pick_rep_set':              {}},
        #{'make_otu_table':            {}},
        #{'filter_otu_table':          {}},
        #{'plot_cluster_distribution': {}},
        #{'plot_otu_distribution':     {}},
        #{'plot_sample_sums':          {}},
        #{'plot_otu_sums':             {}},
        #{'nmds':             {}},
        #{'permanova':             {}},
        {'beta_dispersion':             {}},
    ]

    all_paths = """
    /representatives/rep_set.fasta
    /taxonomy/pool1_hist.pdf
    /taxonomy/pool1_hist.csv
    /taxonomy/pool2_hist.pdf
    /taxonomy/pool2_hist.csv
    /taxonomy/pool3_hist.pdf
    /taxonomy/pool3_hist.csv
    /taxonomy/pool4_hist.pdf
    /taxonomy/pool4_hist.csv
    /taxonomy/pool5_hist.pdf
    /taxonomy/pool5_hist.csv
    /taxonomy/rep_set.taxonomy
    /graphs/cluster_hist.pdf
    /graphs/otu_hist.pdf
    /graphs/sample_sums_hist.pdf
    /graphs/otu_sums_hist.pdf
    /table/table.biom
    /table/table.csv
    /table/table_unfiltered.csv
    /table/table_filtered.csv
    /table/table_trimmed.csv
    /table/table_transposed.csv
    /nmds/by_barcode.pdf
    /nmds/by_pool.pdf
    /nmds/by_run.pdf
    /nmds/r.out
    /permanova/all.txt
    /permanova/pairs.txt
    /beta_dispersion/permutest.txt
    /beta_dispersion/anova.txt
    /beta_dispersion/pool_plots.pdf
    /beta_dispersion/run_plots.pdf
    """

    def __init__(self, base_dir, orig_reads, procedure):
        # Paths #
        if not base_dir.endswith('/') : base_dir += '/'
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Attributes #
        self.orig_reads = orig_reads
        self.procedure = procedure
        self.taxonomy = None

    def pick_rep_set(self):
        pick_rep = sh.Command('pick_rep_set.py')
        pick_rep('-i', self.p.clusters_otus_txt, '-f', self.orig_reads, '-o', self.p.rep_set_fasta)

    def make_otu_table(self):
        # Make BIOM table #
        make_table = sh.Command('make_otu_table.py')
        if self.taxonomy: make_table('-i', self.p.clusters_otus_txt, '-t', self.p.rep_set_taxonomy, '-o', self.p.biom_table)
        else:             make_table('-i', self.p.clusters_otus_txt, '-o', self.p.biom_table)
        # Make CSV table #
        convert_table = sh.Command('convert_biom.py')
        convert_table('-i', self.p.biom_table, '-o', self.p.csv_table, '-b')
        # Remove first line #
        sh.sed('-i', '1d', self.p.csv_table)
        # Remove space in OTU ID #
        sh.sed('-i', '1s/^#OTU ID/OTUID/', self.p.csv_table)

    def filter_otu_table(self):
        # Check sum is at least 3 and back to int #
        def goodlines():
            handle = open(self.p.csv_table)
            yield handle.next()
            for line in handle:
                line = line.split()
                label = line[0]
                values = map(float, line[1:])
                values = map(int, values)
                if sum(values) > 2:
                    yield label + '\t' + '\t'.join(map(str,values)) + '\n'
        handle = open(self.p.csv_filtered_table, 'w')
        handle.writelines(goodlines())
        handle.close()
        # Transpose #
        handle = open(self.p.csv_transposed_table, "w")
        rows = izip(*csv.reader(open(self.p.csv_filtered_table), delimiter='\t'))
        csv.writer(handle, delimiter='\t').writerows(rows)
        handle.close()
        # Remove low samples #
        def goodlines():
            handle = open(self.p.transposed_table)
            yield handle.next()
            for line in handle:
                line = line.split()
                label = line[0]
                values = map(int, line[1:])
                if sum(values) > 10:
                    yield label + '\t' + '\t'.join(map(str,values)) + '\n'
        handle = open(self.p.table_trimmed, 'w')
        handle.writelines(goodlines())
        handle.close()

    def nmds(self):
        # Script #
        script = []
        # Load libs #
        script += ["library(vegan)"]
        script += ["library(MASS)"]
        script += ["library(ggplot2)"]
        script += ["library(compare)"]
        # Load data #
        script += ["data = read.table('%s', header=TRUE, sep='\t', row.names='OTUID')" % (self.p.table_trimmed)]
        script += ["meta = read.table('%s', header=TRUE, sep='\t', row.names=1)" % (self.procedure.metadata_path)]
        # Compute nmds #
        script += ["ord = metaMDS(data, distance='%s', trymax=200)" % self.dist_method]
        script += ["nmds = scores(ord)"]
        # Make a dataframe #
        script += ["df = merge(meta, nmds, by.x='row.names', by.y='row.names')"]
        # Make factors #
        script += ["df$barcode = factor(df$barcode)"]
        script += ["df$run = factor(df$run)"]
        script += ["df$pool = factor(df$pool)"]
        # Make plots #
        script += ["p = ggplot(df, aes(NMDS1, NMDS2)) + xlab('Dimension 1') + ylab('Dimension 2')"]
        script += ["pdf(file='%s')" % self.p.nmds_by_barcode]
        script += ["p + geom_point(aes(colour=barcode))"]
        script += ["dev.off()"]
        script += ["pdf(file='%s')" % self.p.nmds_by_pool]
        script += ["p + geom_point(aes(colour=pool))"]
        script += ["dev.off()"]
        script += ["pdf(file='%s')" % self.p.nmds_by_run]
        script += ["p + geom_point(aes(colour=run))"]
        script += ["dev.off()"]
        # Run it #
        sh.R('--no-save', '-f', TmpFile.from_string('\n'.join(script)), _out=self.p.nmds_out)

    def permanova(self):
        # Basic PERMANOVA #
        ro.r("library(vegan)")
        ro.r("data = read.table('%s', header=TRUE, sep='\t', row.names='OTUID')" % (self.p.table_trimmed))
        ro.r("meta = read.table('%s', header=TRUE, sep='\t', row.names=1)" % (self.procedure.metadata_path))
        ro.r("data_ordered = data[order(row.names(data)),]")
        ro.r("meta_ordered = meta[row.names(data),]")
        ro.r("meta_ordered = meta_ordered[order(row.names(meta_ordered)),]")
        # As factor #
        ro.r("meta_ordered$pool = factor(meta_ordered$pool)")
        ro.r("meta_ordered$barcode = factor(meta_ordered$barcode)")
        ro.r("meta_ordered$run = factor(meta_ordered$run)")
        # Run test #
        ro.r("permanova = adonis(formula = data ~ pool * barcode, data=meta_ordered, permutations=1000, method='%s')" % self.dist_method)
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

    def beta_dispersion(self):
        # Prepare #
        ro.r("library(vegan)")
        ro.r("data = read.table('%s', header=TRUE, sep='\t', row.names='OTUID')" % (self.p.table_trimmed))
        ro.r("meta = read.table('%s', header=TRUE, sep='\t', row.names=1)" % (self.procedure.metadata_path))
        # Compute #
        ro.r("data_ordered = data[order(row.names(data)),]")
        ro.r("data_sqrt = sqrt(data_ordered)")
        ro.r("data_wa = wisconsin(data_sqrt)")
        ro.r("data_dist = vegdist(data_wa, method='%s')" % self.dist_method)
        ro.r("meta = meta[row.names(data),]")
        ro.r("pool = factor(meta[,1])")
        ro.r("run = factor(meta[,3])")
        # Run group #
        ro.r("mod1 = betadisper(data_dist, run)")
        ro.r("test = permutest(mod1, control = permControl(nperm = 1000))")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_permutest, 'w') as handle: handle.write(result)
        ro.r("test = anova(mod1)")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_anova, 'w') as handle: handle.write(result)
        ro.r("pdf(file='%s')" % self.p.beta_dispersion_run_plots)
        ro.r("plot(mod1)")
        ro.r("boxplot(mod1)")
        ro.r("plot(TukeyHSD(mod1))")
        ro.r("dev.off()")
        # Pool group #
        ro.r("mod2 = betadisper(data_dist, pool)")
        ro.r("test = permutest(mod2, control = permControl(nperm = 1000))")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_permutest, 'w') as handle: handle.write(result)
        ro.r("mod1 = betadisper(data_dist, run)")
        ro.r("test = anova(mod2)")
        result = '\n'.join(ro.r("capture.output(print(test))")).encode('utf-8')
        with open(self.p.beta_dispersion_anova, 'w') as handle: handle.write(result)
        ro.r("pdf(file='%s')" % self.p.beta_dispersion_pool_plots)
        ro.r("plot(mod2)")
        ro.r("boxplot(mod2)")
        ro.r("plot(TukeyHSD(mod2))")
        ro.r("dev.off()")

    def plot_cluster_distribution(self): cluster_distribution(self)
    def plot_otu_distribution(self): otu_distribution(self)
    def plot_sample_sums(self): sample_sums(self)
    def plot_otu_sums(self): otu_sums(self)
    def plot_taxa_hist(self): taxa_hist(self)

###############################################################################
class DenovoOTUs(OTUs):
    short_name = 'denovo'
    method = 'Denovo picking'

    all_paths = OTUs.all_paths + """
    /clusters/otus.txt
    /clusters/otus.log
    /clusters/clusters.uc
    """

    def pick_otus(self):
        # Prepare #
        pick_otus = sh.Command('pick_otus.py')
        shutil.rmtree(self.p.clusters_dir)
        # Run command #
        pick_otus('-m', 'uclust', '-s', 0.97, '-i', self.orig_reads, '-o', self.p.clusters_dir)
        # Move into place #
        base_name = self.p.clusters_dir + os.path.basename(self.orig_reads)[:-6]
        shutil.move(base_name + '_clusters.uc', self.p.clusters_uc)
        shutil.move(base_name + '_otus.log', self.p.clusters_otus_log)
        shutil.move(base_name + '_otus.txt', self.p.clusters_otus_txt)

#------------------------------------------------------------------------------#
class OpenRefOTUs(OTUs):
    short_name = 'openref'
    method = 'Open reference picking'

    all_paths = OTUs.all_paths + """
    /clusters/otus.txt
    /clusters/otus.log
    /clusters/clusters.uc
    """

    def pick_otus(self):
        # Prepare #
        pick_otus = sh.Command('pick_open_reference_otus.py')
        shutil.rmtree(self.p.otus_dir)
        green_genes_db_path = home + "share/green_genes/rep_set/97_otus.fasta"
        # Run command #
        pick_otus('-m', 'uclust', '-i', self.orig_reads, '-o', self.p.clusters_dir, '-f', '-a', '-O', 8,
                  '--reference_fp', green_genes_db_path,
                  '-p', TmpFile.from_string('pick_otus:enable_rev_strand_match False'),
                  '--suppress_align_and_tree')

#------------------------------------------------------------------------------#
class StepOTUs(OTUs):
    short_name = 'stepwise'
    method = '99-98-97 progre   ssive picking'

    all_paths = OTUs.all_paths + """
    /clusters/otus_99.txt
    /clusters/otus_99.log
    /clusters/clusters_99.uc
    /clusters/rep_set_99.fasta
    /clusters/otus_98.txt
    /clusters/otus_98.log
    /clusters/clusters_98.uc
    /clusters/rep_set_98.fasta
    /clusters/otus_97.txt
    /clusters/otus_97.log
    /clusters/clusters_97.uc
    /clusters/rep_set_97.fasta
    /clusters/otus.txt
    """

    def pick_otus(self):
        # Commands #
        pick_otus = sh.Command('pick_otus.py')
        pick_rep = sh.Command('pick_rep_set.py')
        # 99 #
        pick_otus('-m', 'uclust', '-s', 0.99, '-i', self.orig_reads, '-o', self.p.clusters_dir)
        base_name = self.p.clusters_dir + os.path.basename(self.orig_reads)[:-6]
        shutil.move(base_name + '_clusters.uc', self.p.clusters_99_uc)
        shutil.move(base_name + '_otus.log', self.p.otus_99_log)
        shutil.move(base_name + '_otus.txt', self.p.otus_99_txt)
        pick_rep('-i', self.p.otus_99_txt, '-f', self.orig_reads, '-o', self.p.rep_set_99_fasta)
        # 98 #
        pick_otus('-m', 'uclust', '-s', 0.98, '-i', self.p.rep_set_99_fasta, '-o', self.p.clusters_dir)
        base_name = self.p.clusters_dir + os.path.basename(self.p.rep_set_99_fasta)[:-6]
        shutil.move(base_name + '_clusters.uc', self.p.clusters_98_uc)
        shutil.move(base_name + '_otus.log', self.p.otus_98_log)
        shutil.move(base_name + '_otus.txt', self.p.otus_98_txt)
        pick_rep('-i', self.p.otus_98_txt, '-f', self.p.rep_set_99_fasta, '-o', self.p.rep_set_98_fasta)
        # 97 #
        pick_otus('-m', 'uclust', '-s', 0.97, '-i', self.p.rep_set_98_fasta, '-o', self.p.clusters_dir)
        base_name = self.p.clusters_dir + os.path.basename(self.p.rep_set_98_fasta)[:-6]
        shutil.move(base_name + '_clusters.uc', self.p.clusters_97_uc)
        shutil.move(base_name + '_otus.log', self.p.otus_97_log)
        shutil.move(base_name + '_otus.txt', self.p.otus_97_txt)
        pick_rep('-i', self.p.otus_97_txt, '-f', self.p.rep_set_98_fasta, '-o', self.p.rep_set_97_fasta)
        # Read children #
        childs_99 = {}
        for line in open(self.p.otus_99_txt):
            line = line.strip('\n').split()
            childs_99[line.pop(0)] = line
        childs_98 = {}
        for line in open(self.p.otus_98_txt):
            line = line.strip('\n').split()
            childs_98[line.pop(0)] = line
        childs_97 = {}
        for line in open(self.p.otus_97_txt):
            line = line.strip('\n').split()
            childs_97[line.pop(0)] = line
        # Combine children #
        clusters = {}
        for key in childs_97:
            reads = flatten([childs_98[v] for v in childs_97[key]])
            reads = flatten([childs_99[v] for v in reads])
            clusters[key] = reads
        # Combine children #
        with open(self.p.otus_txt, 'w') as handle:
            for k,v in clusters.items(): handle.write(k + '\t' + '\t'.join(v) + '\n')

#------------------------------------------------------------------------------#
class SubsampledOTUs(OTUs):
    short_name = 'subsampled'
    method = 'Denovo Subsampled OTUs'
    dist_method = 'bray'

    all_paths = OTUs.all_paths + """
    /table/subsampled_float.csv
    """

    default_steps = [
        #{'subsample_table':           {}},
        #{'plot_sample_sums':          {}},
        #{'plot_otu_sums':             {}},
        #{'nmds':                      {}},
        #{'permanova':                      {}},
        {'beta_dispersion':                      {}},
    ]

    def __init__(self, base_dir, base_otu, procedure):
        # Paths #
        if not base_dir.endswith('/') : base_dir += '/'
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Attributes #
        self.base_otu = base_otu
        self.procedure = procedure

    def filter_otu_table(self): raise NotImplementedError('')

    def subsample_table(self):
        # Parse #
        otus = pandas.read_csv(self.base_otu.p.table_trimmed, sep = '\t', index_col=0)
        # Drop those below 10 #
        #sums = otus.sum(axis=1)
        #otus = otus.drop(sums[sums < 10].keys())
        # Subsample #
        sums = otus.sum(axis=1)
        down_to = min(sums)
        subotus = pandas.DataFrame(columns=otus.columns, index=otus.index, dtype=int)
        # Do it #
        for sample_name in otus.index:
            row = otus.loc[sample_name]
            weighted_choices = list(row[row != 0].iteritems())
            population = [val for val, count in weighted_choices for i in range(count)]
            sub_pop = random.sample(population, down_to)
            frequencies = Counter(sub_pop)
            new_row = pandas.Series(frequencies.values(), index=frequencies.keys(), dtype=int)
            subotus.loc[sample_name] = new_row
        # Output it #
        subotus.to_csv(self.p.subsampled_table_float, sep='\t', na_rep='0')
        # Cast to integer #
        def lines_as_integer(path):
            handle = open(path)
            yield handle.next()
            for line in handle:
                line = line.split()
                label = line[0]
                values = map(float, line[1:])
                values = map(int, values)
                yield label + '\t' + '\t'.join(map(str,values)) + '\n'
        handle = open(self.p.table_trimmed, 'w')
        handle.writelines(lines_as_integer(self.p.subsampled_table_float))
        handle.close()