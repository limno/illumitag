# Built-in modules #
import os

# Internal modules #
from illumitag.common import flatten
from illumitag.common.autopaths import AutoPaths
from illumitag.common.tmpstuff import TmpFile
from illumitag.common.csv_tables import TSVTable
from illumitag.analysis.statistics import StatsOnOTU
from illumitag.graphs import otu_plots

# Third party modules #
import shutil, sh

# Constants #
home = os.environ['HOME'] + '/'
green_genes_db_path = home + "share/green_genes/rep_set/97_otus.fasta"

###############################################################################
class OTUs(object):
    """Base class for OTU analyses"""

    all_paths = """
    /representatives/rep_set.fasta
    /table/table.biom
    /table/table.csv
    /table/table_filtered.csv
    /table/table_transposed.csv
    /graphs/
    /subsampled/
    """

    def __repr__(self): return '<%s object %s of %s>' % \
                               (self.__class__.__name__, self.short_name, self.parent)

    def __init__(self, parent):
        # Save parent #
        self.analysis, self.parent = parent, parent
        # Inherited #
        self.pools = self.parent.pools
        self.qiime_reads = self.parent.qiime_reads
        self.meta_data_path = self.parent.meta_data_path
        # Paths #
        self.base_dir = self.parent.p.otus_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Other #
        self.taxonomy = None
        # Files #
        self.table = TSVTable(self.p.csv_table)
        self.table_filtered = TSVTable(self.p.csv_table_filtered)
        self.table_transposed = TSVTable(self.p.csv_table_transposed)
        # Children #
        self.stats = StatsOnOTU(self, self.table_filtered)
        # Deferred import #
        from illumitag.analysis.subsample import SubsampledOTUs
        self.subsampled = SubsampledOTUs(self)

    def run(self):
        # Standard #
        self.pick_otus()
        self.pick_rep_set()
        self.make_otu_table()
        self.filter_otu_table()
        self.make_otu_plots()
        self.compute_stats()
        # Subsample #
        self.subsampled.run()

    def pick_rep_set(self):
        pick_rep = sh.Command('pick_rep_set.py')
        pick_rep('-i', self.p.clusters_otus_txt, '-f', self.qiime_reads.path, '-o', self.p.rep_set_fasta)

    def make_otu_table(self):
        # Make BIOM table #
        make_table = sh.Command('make_otu_table.py')
        if self.taxonomy: make_table('-i', self.p.clusters_otus_txt, '-t', self.p.rep_set_taxonomy, '-o', self.p.biom_table)
        else:             make_table('-i', self.p.clusters_otus_txt, '-o', self.p.biom_table)
        # Make CSV table #
        convert_table = sh.Command('convert_biom.py')
        convert_table('-i', self.p.biom_table, '-o', self.table.path, '-b')

    def filter_otu_table(self):
        # Format it #
        self.table.remove_first_line()
        self.table.replace_title('#OTU ID', 'OTUID')
        self.table.to_integer()
        # Make a filtered version #
        shutil.copy(self.table.path, self.table_filtered.path)
        self.table_filtered.filter_line_sum(minimum=3) # Min cluster size
        self.table_filtered.transpose()
        self.table_filtered.filter_line_sum(minimum=10) # Min reads in sample
        # Make a transposed version #
        self.table_filtered.transpose(path=self.table_transposed.path)

    def compute_stats(self):
        self.stats.run()

    def make_otu_plots(self):
        for cls_name in otu_plots.__all__:
            cls = getattr(otu_plots, cls_name)
            cls(self).plot()

###############################################################################
class DenovoOTUs(OTUs):
    short_name = 'denovo'
    method = 'Denovo picking'
    dist_method = 'horn'

    all_paths = OTUs.all_paths + """
    /clusters/otus.txt
    /clusters/otus.log
    /clusters/clusters.uc
    """

    def pick_otus(self):
        # Clean #
        shutil.rmtree(self.p.clusters_dir)
        # Prepare #
        pick_otus = sh.Command('pick_otus.py')
        # Run command #
        pick_otus('-m', 'uclust', '-s', 0.97, '-i', self.qiime_reads.path, '-o', self.p.clusters_dir)
        # Move into place #
        base_name = self.p.clusters_dir + self.qiime_reads.prefix
        shutil.move(base_name + '_otus.txt', self.p.clusters_otus_txt)
        shutil.move(base_name + '_otus.log', self.p.clusters_otus_log)
        shutil.move(base_name + '_clusters.uc', self.p.clusters_uc)

###############################################################################
class OpenRefOTUs(OTUs):
    short_name = 'openref'
    method = 'Open reference picking'
    dist_method = 'horn'

    all_paths = OTUs.all_paths + """
    /clusters/otus.txt
    /clusters/otus.log
    /clusters/clusters.uc
    """

    def pick_otus(self):
        # Clean #
        shutil.rmtree(self.p.clusters_dir)
        # Prepare #
        pick_otus = sh.Command('pick_open_reference_otus.py')
        # Run command #
        pick_otus('-m', 'uclust', '-i', self.orig_reads, '-o', self.p.clusters_dir, '-f', '-a', '-O', 8,
                  '--reference_fp', green_genes_db_path,
                  '-p', TmpFile.from_string('pick_otus:enable_rev_strand_match False'),
                  '--suppress_align_and_tree')
        # Move into place #
        base_name = self.p.clusters_dir + self.qiime_reads.prefix
        shutil.move(base_name + '_otus.txt', self.p.clusters_otus_txt)
        shutil.move(base_name + '_otus.log', self.p.clusters_otus_log)
        shutil.move(base_name + '_clusters.uc', self.p.clusters_uc)

###############################################################################
class StepOTUs(OTUs):
    short_name = 'stepwise'
    method = '99-98-97 progressive picking'
    dist_method = 'horn'

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
        # Clean #
        shutil.rmtree(self.p.clusters_dir)
        # Commands #
        pick_otus = sh.Command('pick_otus.py')
        pick_rep = sh.Command('pick_rep_set.py')
        # 99 #
        pick_otus('-m', 'uclust', '-s', 0.99, '-i', self.qiime_reads.path, '-o', self.p.clusters_dir)
        base_name = self.p.clusters_dir + self.qiime_reads.prefix
        shutil.move(base_name + '_clusters.uc', self.p.clusters_99_uc)
        shutil.move(base_name + '_otus.log', self.p.otus_99_log)
        shutil.move(base_name + '_otus.txt', self.p.otus_99_txt)
        pick_rep('-i', self.p.otus_99_txt, '-f', self.qiime_reads.path, '-o', self.p.rep_set_99_fasta)
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