# Built-in modules #
import os, time, shutil, tempfile

# Internal modules #
from illumitag.cluster.taxonomy import Taxonomy
from illumitag.common.cache import property_cached
from illumitag.helper.mothur import process_log_file

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #
train_db_path = home + 'glob/16s/trainset/trainset9_032012.pds.fasta'
train_tax_path = home + 'glob/16s/trainset/trainset9_032012.pds.tax'
silva_db_path = home + "/share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta"

###############################################################################
class MothurTaxonomy(Taxonomy):
    def assign(self):
        # Prepare #
        self.tmp_dir = tempfile.mkdtemp() + '/'
        os.symlink(self.fasta.path, self.tmp_dir + 'reads.fasta')
        # Run #
        start_time = time.asctime()
        sh.mothur("#classify.seqs(fasta=%s, reference=%s, taxonomy=%s, cutoff=%i, processors=8, probs=F)" %
                 (self.tmp_dir + 'reads.fasta', train_db_path, train_tax_path, 80))
        process_log_file('mothur.classify.seqs.logfile', self.base_dir, start_time)

    @property_cached
    def assignments(self):
        # Save the flipped ids #
        flipped_ids = [line.strip('\n') for line in open(self.tmp_dir + 'reads.pds.wang.flip.accnos')]
        # Add sample name via hash table #
        id_to_sample_hash = dict((r.description.split() for r in self.fasta.parse()))
        # Remove stupid quotes #
        with open(self.p.reads_taxonomy, 'w') as handle:
            for line in open(self.tmp_dir + 'reads.pds.wang.taxonomy', 'r'):
                the_id, taxon = line.split()
                if the_id in flipped_ids: taxon = "unknown;unclassified;unclassified;unclassified;unclassified;unclassified;"
                handle.write(the_id + ' ' + id_to_sample_hash[the_id] + '\t' + taxon.replace('"','') + '\n')
        # Cleanup #
        shutil.rmtree(self.tmp_dir)
