# Built-in modules #
import os, time, shutil, tempfile

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA
from illumitag.helper.mothur import process_log_file
from illumitag.common.slurm import nr_threads
from illumitag.common.cache import property_cached

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

# Databases #
train_db_path = home + 'glob/16s/trainset/trainset9_032012.pds.fasta'
train_tax_path = home + 'glob/16s/trainset/trainset9_032012.pds.tax'
silva_db_path = home + "/share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta"

###############################################################################
class Taxonomy(object):
    """Can assign taxonomy to a FASTA file of 16S sequences."""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, fasta_path, parent):
        # Parent #
        self.otu, self.parent = parent, parent
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # Dir #
        self.base_dir = self.parent.p.taxonomy_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def at_level(self, level):
        return dict((k,v[level]) for k,v in self.assignments)

###############################################################################
class CrestTaxonomy(Taxonomy):
    all_paths = """
    /reads.taxonomy
    /silva_hits.xml
    /silva_composition.txt
    /silva_tree.txt
    /silva_assignments.txt
    """

    def assign(self):
        # Run #
        sh.megablast('-a', nr_threads, '-i', self.fasta, '-d', silva_db_path, '-b100', '-v100', '-m7', '-o', self.p.silva_hits)
        if os.path.getsize(self.p.silva_hits) == 0: raise Exception("Hits file empty. The MEGABLAST process was probably killed.")
        # CREST #
        sh.classify(self.p.silva_hits, '-p', '-o', '-d', 'silvamod')
        shutil.move(self.p.silva_hits[:-4] + '_Composition.txt', self.p.silva_composition)
        shutil.move(self.p.silva_hits[:-4] + '_Tree.txt', self.p.silva_tree)
        shutil.move(self.p.silva_hits[:-4] + '_Assignments.txt', self.p.silva_assignments)

    @property_cached
    def assignments(self):
        result = {}
        for line in open(self.p.assignments):
            code, species = line.split('\t')
            result[code] = tuple(species.strip('\n').split(';'))[:8]
        return result

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