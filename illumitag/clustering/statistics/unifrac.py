# Built-in modules #
import os

# Internal modules #
from illumitag.common.autopaths import AutoPaths

# Third party modules #
import cogent, sh
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_tree import UniFracTreeNode
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
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
        * mothur
        * SINA

    Step 2. From the alignment make a tree
    Step 3. Feed the tree and the otu table to pycogent
    Step 4. Return distance matrix

    http://telliott99.blogspot.se/2010/02/unifrac-analysis-introduction.html"""

    all_paths = """
    /clustalo/centers.align
    /pynast/centers.align
    /pynast/log.txt
    /pynast/fail.fasta
    /mothur/
    /sina/
    /centers.align
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

    def run(self):
        # Step 1 #
        self.align_mothur()
        # Step 2 #
        pass
        # Step 3 #
        pass
        # Step 4 #
        pass

    def align_clustalo(self):
        # Step 1 clustalo #
        self.centers_aligned.remove()
        sh.clustalo('-i', self.tax.centers, '--profile1', reference, '-o', self.clustalo_aligned, '--threads', 16)

    def align_pynast(self):
        # Step 1 PyNAST #
        sh.pynast('--input_fp', self.tax.centers, '--template_fp', reference,
                  '--fasta_out_fp', self.pynast_aligned,
                  '--log_fpa', self.p.pynast_log, '--failure_fp', self.p.pynast_fail,)

    def align_mothur(self):
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=blast, flip=false, processors=16);" % (self.path, ref_path))
        # Move things #
        shutil.move(self.path[:-6] + '.align', self.aligned_path)
        shutil.move(self.path[:-6] + '.align.report', self.report_path)
        shutil.move(self.path[:-6] + '.flip.accnos', self.accnos_path)
        # Clean up #
        if os.path.exists('formatdb.log'): os.remove('formatdb.log')
        if os.path.exists('error.log') and os.path.getsize('error.log') == 0: os.remove('error.log')
        for p in sh.glob('mothur.*.logfile'): os.remove(p)

    def align_sina(self):
        # Step 1 SINA: should use the NR 99% silva database and PT-SERVER #
        pass


