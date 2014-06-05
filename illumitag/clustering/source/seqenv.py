# Built-in modules #
import os, shutil

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.common.slurm import nr_threads
from illumitag.common.csv_tables import CSVTable
from illumitag.common.tmpstuff import TmpFile

# Third party modules #
import pandas, sh
from Bio import SeqIO

# Constants #
home = os.environ['HOME'] + '/'
seqenv_script = home + "share/seqenv/SEQenv_v0.8/SEQenv_samples.sh"

###############################################################################
class Seqenv(object):
    """Base class for Seqenv results processing."""
    N = 1000

    all_paths = """
    /working_dir/
    /abundances.csv
    /seqenv.out
    /
    /centers_N1000_blast_F_ENVO_OTUs_labels.csv
    """

    def __init__(self, parent, base_dir=None):
        # Parent #
        self.otu, self.parent = parent, parent
        self.taxonomy = self.parent.taxonomy
        # Inherited #
        self.samples = self.parent.samples
        # Dir #
        if base_dir is None: self.base_dir = self.parent.p.seqenv
        else: self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        self.abundances = CSVTable(self.p.abundances)

    def tr(self): self.taxonomy.otu_csv_norm.transpose(self.abundances, d=',')

    def run(self):
        # Move to the working dir #
        self.saved_cwd = os.getcwd()
        os.chdir(self.p.working_dir)
        # Make the abundances file #
        self.taxonomy.otu_csv_norm.transpose(self.abundances, d=',')
        # Make the most abundant OTU file (launching one perl command per OTU sequence takes forever) #
        path = "centers_N%i.fa" % self.N
        if len(self.taxonomy.centers) <= self.N: self.taxonomy.centers.copy(path)
        else: # Untested
            otus = self.taxonomy.otu_table_norm.sum()
            otus.sort()
            highest_otus = otus[0:self.N]
            sequences = (seq for seq in SeqIO.parse(self.taxonomy.centers, 'fasta') if seq.id in highest_otus)
            with open(path, 'w') as handle: SeqIO.write(sequences, handle, 'fasta')
        # Run the Quince pipeline with a special version of R #
        module = "module load R/3.0.1"
        params = ['-f', self.taxonomy.centers, '-s', self.abundances, '-n', self.N, '-p', '-c', nr_threads]
        command = "bash -x " + seqenv_script + ' ' + ' '.join(map(str,params))
        sh_file = TmpFile.from_string(module + '\n' + command)
        sh.bash(sh_file, _out=str(self.p.out))
        # Move things into place #
        shutil.move("centers_N1000_blast_F_ENVO_OTUs.csv", "../")
        shutil.move("centers_N1000_blast_F_ENVO_OTUs_labels.csv", "../")
        # Cleanup #
        os.chdir(self.saved_cwd)
        #shutil.rmtree(self.p.working_dir)

    @property
    def frame(self):
        return pandas.io.parsers.read_csv(self.p.labels, sep=',', index_col=0, encoding='utf-8')