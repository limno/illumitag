# Built-in modules #
import os, csv
from itertools import izip

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.common.tmpstuff import TmpFile
from illumitag.common import flatten

# Third party modules #
import shutil, sh

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class OTUs(object):
    """Base class for OTU analyses"""

    all_paths = """
    /representatives/rep_set.fasta
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
    """

    def __init__(self, parent):
        # Save parent #
        self.analysis, self.parent = parent, parent
        self.pools = parent.pools
        # Paths #
        self.base_dir = self.parent.p.otus_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

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
        # Prepare #
        pick_otus = sh.Command('pick_open_reference_otus.py')
        shutil.rmtree(self.p.otus_dir)
        green_genes_db_path = home + "share/green_genes/rep_set/97_otus.fasta"
        # Run command #
        pick_otus('-m', 'uclust', '-i', self.orig_reads, '-o', self.p.clusters_dir, '-f', '-a', '-O', 8,
                  '--reference_fp', green_genes_db_path,
                  '-p', TmpFile.from_string('pick_otus:enable_rev_strand_match False'),
                  '--suppress_align_and_tree')

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