# Built-in modules #
import os, time, shutil, tempfile

# Internal modules #
from illumitag.common import AutoPaths
from illumitag.fasta.single import FASTA

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'
train_db_path = home + 'glob/16s/trainset/trainset9_032012.pds.fasta'
train_tax_path = home + 'glob/16s/trainset/trainset9_032012.pds.tax'
silva_path = home + "/share/LCAClassifier/parts/flatdb/silvamod/silvamod.fasta"

###############################################################################
class Classifier(object):
    """Can assign taxonomy to a FASTA file."""

    all_paths = """
    /reads.taxonomy
    /silva_hits.xml
    /silva_composition.txt
    /silva_tree.txt
    /silva_assignments.txt
    """

    def __repr__(self): return '<%s object on "%s">' % (self.fasta.path)

    def __init__(self, fasta_path, base_dir):
        # FASTA #
        self.fasta = FASTA(fasta_path)
        # Dir #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def assign_mothur(self):
        # Prepare #
        directory = tempfile.mkdtemp() + '/'
        os.symlink(self.fasta.path, directory + 'reads.fasta')
        # Run #
        start_time = time.asctime()
        sh.mothur("#classify.seqs(fasta=%s, reference=%s, taxonomy=%s, cutoff=%i, processors=8, probs=F)" %
                 (directory + 'reads.fasta', train_db_path, train_tax_path, 80))
        process_log_file('mothur.classify.seqs.logfile', self.base_dir, start_time)
        # Save the flipped ids #
        flipped_ids = [line.strip('\n') for line in open(directory + 'reads.pds.wang.flip.accnos')]
        # Add sample name via hash table #
        otu_to_sample_hash = dict((r.description.split() for r in self.fasta.parse()))
        # Remove stupid quotes #
        with open(self.p.rep_set_taxonomy, 'w') as handle:
            for line in open(directory + 'reads.pds.wang.taxonomy', 'r'):
                name, taxon = line.split()
                if name in flipped_ids: taxon = "unknown;unclassified;unclassified;unclassified;unclassified;unclassified;"
                handle.write(name + ' ' + otu_to_sample_hash[name] + '\t' + taxon.replace('"','') + '\n')
        # Cleanup #
        shutil.rmtree(directory)

    def assign_crest(self):
        # Run #
        sh.megablast('-a', 8, '-i', self.p.rep_set_fasta, '-d', silva_path, '-b100', '-v100', '-m7', '-o', self.p.silva_hits)
        if os.path.getsize(self.p.silva_hits) == 0: raise Exception("The MEGABLAST process was probably killed. Hits file empty.")
        # CREST #
        sh.classify(self.p.silva_hits, '-p', '-o', '-d', 'silvamod')
        shutil.move(self.p.silva_hits[:-4] + '_Composition.txt', self.p.silva_composition)
        shutil.move(self.p.silva_hits[:-4] + '_Tree.txt', self.p.silva_tree)
        shutil.move(self.p.silva_hits[:-4] + '_Assignments.txt', self.p.silva_assignments)
        # Assignments #
        self.silva_assignments = {}
        for line in open(self.p.silva_assignments):
            code, species = line.split('\t')
            self.silva_assignments[code] = tuple(species.strip('\n').split(';'))[:8]

    def assign_qiime(self):
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
