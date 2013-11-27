# Built-in modules #
import shutil

# Internal modules #
from illumitag.cluster.taxonomy import Taxonomy
from illumitag.common.cache import property_cached

# Third party modules #
import sh

###############################################################################
class QiimeTaxonomy(Taxonomy):
    def assign(self):
        shutil.rmtree(self.p.taxonomy_dir)
        # Empty output #
        assign = sh.Command('assign_taxonomy.py')
        assign('-i', self.p.rep_set_fasta, '-o', self.p.taxonomy_dir, '--rdp_max_memory', 18000)
        # Hangs #
        assign = sh.Command('parallel_assign_taxonomy_rdp.py')
        assign('-i', self.p.rep_set_fasta, '-o', self.p.taxonomy_dir, '--rdp_max_memory', 2500)
        # Crashes #
        assign = sh.Command('assign_taxonomy.py')
        assign('-i', self.p.rep_set_fasta, '-o', self.p.taxonomy_dir, '-m', 'mothur')
        # Missing files #
        assign = sh.Command('assign_taxonomy.py')
        assign('-i', self.p.rep_set_fasta, '-o', self.p.taxonomy_dir, '-m', 'rtax',
               '--read_1_seqs_fp', self.p.p, '--read_2_seqs_fp', self.p.p)
        # Too slow #
        assign = sh.Command('assign_taxonomy.py')
        assign('-i', self.p.rep_set_fasta, '-o', self.p.taxonomy_dir, '-m', 'blast')

    @property_cached
    def assignments(self):
        pass