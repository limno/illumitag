# Built-in modules #
import os, shutil, glob

# Internal modules #
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.common.cache import property_cached
from illumitag.common.csv_tables import CSVTable
from illumitag.clustering.statistics.nmds import NMDS

# Third party modules #
import sh, pandas
from illumitag.fasta.single import FASTA

# Constants #
home = os.environ['HOME'] + '/'
reference = FASTA(home + 'glob/16s/silva/v111/rep_set_aligned/97_Silva_111_rep_set.fasta')

###############################################################################
class Unifrac(object):
    """A class to compute the Unifrac algorithm producing a distance matrix
    from a bunch of different samples and their reads.

    Step 1. Make an alignment of all the OTU centers against a reference.
    One can use:
        * clustalo
        * PyNAST
        * mothur <- fastest
        * SINA

    Step 2. From the alignment produced make a tree.
    One can use:
        * RAxML
        * Fasttree <- fastest

    Step 3. Feed the tree and the otu table into a unifrac algorithm
    One can use:
        * Pycogent

    Step 4. Return distance matrix

    http://telliott99.blogspot.se/2010/02/unifrac-analysis-introduction.html"""

    all_paths = """
    /clustalo/centers.align
    /pynast/centers.align
    /pynast/log.txt
    /pynast/fail.fasta
    /mothur/centers.align
    /mothur/report.txt
    /mothur/flip.accnos
    /raxml/output.tree
    /fasttree/output.tree
    /distances.csv
    /nmds/
    """

    def __init__(self, parent):
        # Save parent #
        self.stat, self.parent = parent, parent
        self.tax = parent.tax
        # Paths #
        self.p = AutoPaths(self.parent.p.unifrac_dir, self.all_paths)
        # Files #
        self.clustalo_aligned = FASTA(self.p.clustalo_align)
        self.pynast_aligned = FASTA(self.p.pynast_align)
        self.mothur_aligned = FASTA(self.p.mothur_align)
        self.raxml_tree = FilePath(self.p.raxml_tree)
        self.fasttree_tree = FilePath(self.p.fasttree_tree)
        self.distances_csv = CSVTable(self.p.distances_csv)
        # Graphs #
        self.nmds = NMDS(self, self.distances_csv, calc_distance=False)

    def run(self):
        # Step 1 #
        self.align_mothur()
        # Step 2 #
        self.tree_fasttree()
        # Step 3 #
        self.unifrac_pycogent()
        # Step 4 #
        return self.distances

    def align_clustalo(self):
        """Step 1 with clustalo"""
        self.centers_aligned.remove()
        sh.clustalo('-i', self.tax.centers, '--profile1', reference, '-o', self.clustalo_aligned, '--threads', 16)

    def align_pynast(self):
        """Step 1 with PyNAST"""
        sh.pynast('--input_fp', self.tax.centers, '--template_fp', reference,
                  '--fasta_out_fp', self.pynast_aligned,
                  '--log_fpa', self.p.pynast_log, '--failure_fp', self.p.pynast_fail,)

    def align_mothur(self):
        """Step 1 with mothur"""
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=kmer, flip=false, processors=%s);" \
                  % (self.tax.centers, reference, 16))
        # Move things #
        shutil.move(self.tax.centers.prefix_path + '.align', self.mothur_aligned)
        shutil.move(self.tax.centers.prefix_path + '.align.report', self.p.mothur_report)
        path = self.tax.centers.prefix_path + '.flip.accnos'
        if os.path.exists(path): shutil.move(path, self.p.mothur_accnos)
        # Clean up #
        for p in glob.glob('mothur.*.logfile'): os.remove(p)

    def align_sina(self):
        """Step 1 with SINA: should use the NR 99% silva database and PT-SERVER ?"""
        sh.SINA('-h')

    def tree_raxml(self):
        """Step 2 with RAxML"""
        sh.raxml('-T', 16, '-s', self.mothur_aligned, '-n', self.raxml_tree, '-m', 'LOREM')

    def tree_fasttree(self):
        """Step 2 with FastTree"""
        sh.FastTreeMP('-fastest' ,'-out', self.fasttree_tree, '-nt', self.mothur_aligned)

    def unifrac_pycogent(self):
        """Step 3 with Pycogent"""
        tree_newick = open(self.fasttree_tree, 'r').read()
        from cogent.parse.tree import DndParser
        from cogent.maths.unifrac.fast_tree import UniFracTreeNode
        tree = DndParser(tree_newick, UniFracTreeNode)
        from cogent.maths.unifrac.fast_unifrac import fast_unifrac
        distances = fast_unifrac(tree, self.tax.otu_table.to_dict())
        # Make a dataframe #
        names = distances['distance_matrix'][1]
        df = pandas.DataFrame(distances['distance_matrix'][0], index=names, columns=names)
        df.to_csv(self.distances_csv, sep='\t', float_format='%.5g')

    @property_cached
    def distances(self):
        return pandas.io.parsers.read_csv(self.distances_csv, sep='\t', index_col=0)
