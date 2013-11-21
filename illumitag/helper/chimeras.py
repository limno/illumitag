# Futures #
from __future__ import division

# Built-in modules #
import sys, os, time

# Internal modules #
from illumitag.fasta.single import FASTA, SizesFASTA
from illumitag.common.autopaths import AutoPaths
from illumitag.common.color import Color

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'
chimera_ref_path = home + 'glob/16s/microbiomeutil-r20110519.fasta'

################################################################################
class ChimerasChecker(object):

    all_paths = """
    /subsampled.fasta
    /cluster_99.fasta
    /derep_cluster.fasta
    /positive.fasta
    /negative.fasta
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.fasta)
    def __str__(self): return "%s on %s of pool %s" % (self.long_name, self.parent.parent.parent.long_name,
                                                       self.parent.parent.parent.parent.id_name)

    def __init__(self, fasta_path, base_dir, parent, verbose=False):
        # Base #
        self.fasta = FASTA(fasta_path)
        self.parent = parent
        self.verbose = verbose
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        self.derep_cluster = SizesFASTA(self.p.derep_cluster)
        self.cluster_99 = SizesFASTA(self.p.cluster_99)
        self.positive = SizesFASTA(self.p.positive)
        self.negative = SizesFASTA(self.p.negative)
        self.subsampled = FASTA(self.p.subsampled)

    def clean(self):
        # Clean output #
        open(self.derep_cluster.path, 'w').close()
        open(self.cluster_99.path, 'w').close()
        open(self.positive.path, 'w').close()
        open(self.negative.path, 'w').close()
        open(self.subsampled.path, 'w').close()

    @property
    def percent(self):
        if not self.fasta: return -1
        return (len(self.positive) / len(self.subsampled)) * 100.0

################################################################################
class UchimeRef(ChimerasChecker):
    downto = 100000
    title = "UCHIME algorithm for downsampled pool %i (in reference mode)"
    long_name = "UCHIME reference"

    def check(self):
        # Prepare #
        self.clean()
        if self.verbose: print Color.l_ylw + "----> Subsampling down to %i" % (self.downto,) + Color.end
        sys.stdout.flush()
        self.fasta.subsample(self.downto, self.subsampled.path)
        # Cluster #
        if self.verbose: print Color.l_ylw + "----> Running derep_prefix on %s" % (self.subsampled.path,) + Color.end
        sys.stdout.flush()
        command = ("-derep_prefix", self.subsampled.path, "-output", self.derep_cluster.path, "-sizeout")
        sh.usearch(command)
        # Detect #
        if self.verbose: print Color.l_ylw + "----> Running uchime_ref on %s" % (self.derep_cluster.path,) + Color.end
        sys.stdout.flush()
        command = ("-uchime_ref", self.derep_cluster.path, "-db", chimera_ref_path, "-strand", "plus",
                   "-chimeras", self.positive.path, "-nonchimeras", self.negative.path)
        sh.usearch(command)

################################################################################
class UchimeDenovo(ChimerasChecker):
    downto = 50000
    title = "UCHIME algorithm for downsampled pool %i (in denovo mode)"
    long_name = "UCHIME denovo"

    def check(self):
        # Prepare #
        self.clean()
        if self.verbose: print Color.l_ylw + "----> Subsampling down to %i" % (self.downto,) + Color.end
        sys.stdout.flush()
        self.fasta.subsample(self.downto, self.subsampled.path)
        # Cluster #
        if self.verbose: print Color.l_ylw + "----> Running cluster_fast on %s" % (self.subsampled.path,) + Color.end
        sys.stdout.flush()
        command = ("-cluster_fast", self.subsampled.path, "-id", '0.99',
                   "-consout", self.cluster_99.path, "-sizeout")
        sh.usearch(command)
        # Detect #
        if self.verbose: print Color.l_ylw + "----> Running uchime_denovo on %s" % (self.cluster_99.path,) + Color.end
        sys.stdout.flush()
        command = ("-uchime_denovo", self.cluster_99.path,
                   "-chimeras", self.positive.path, "-nonchimeras", self.negative.path)
        sh.usearch(command)

################################################################################
class Mothur(ChimerasChecker):

    all_paths = ChimerasChecker.all_paths + """
    /asdfasdfasdf
    """

    def __init__(self, *args, **kwargs):
        # Super #
        ChimerasChecker.__init__(self, *args, **kwargs)
        # Other #
        pass

    def check(self):
        start_time = time.asctime()
        sh.mothur("#chimera.uchime(fasta=%s, dereplicate=t, processors=8)" %
                 (self.p.len_filtered_fasta))
        process_log_file('mothur.chimera.uchime.logfile', self.p.orig_dir, start_time)
        #reads.uchime.pick.count_table
        #reads.uchime.chimeras
        #reads.uchime.accnos

################################################################################
class ChimerasSlayer(ChimerasChecker):
    short_name = "chimeras_slayer"
    title = "ChimerasSlayer algorithm for downsampled pool %i"
    downto = 10

    all_paths = ChimerasChecker.all_paths + """
    /subsampled.fasta.NAST.CPS
    /subsampled.fasta.NAST.CPS.CPC
    /subsampled.fasta.NAST.CPS.CPC.wTaxons
    """

    def check(self):
        # Prepare #
        self.clean()
        print Color.l_ylw + "----> Subsampling down to %i" % (self.downto,) + Color.end
        sys.stdout.flush()
        self.fasta.subsample(self.downto, self.subsampled.path)
        # Detect #
        print Color.l_ylw + "----> Running slayer on %s" % (self.subsampled.path,) + Color.end
        sys.stdout.flush()
        slayer = sh.Command(home + "share/microbiomeutil/ChimeraSlayer/ChimeraSlayer.pl")
        slayer('--query_FASTA', self.subsampled.path)

################################################################################
class Perseus(ChimerasChecker):

    all_paths = ChimerasChecker.all_paths + """
    /asdfasdfasd
    """

    def check(self, downto=2):
        # Prepare #
        self.clean()
        print Color.l_ylw + "----> Subsampling down to %i" % (downto,) + Color.end
        sys.stdout.flush()
        self.fasta.subsample(downto, self.subsampled.path)
        # Detect #
        sh.Perseus('-sin', self.subsampled.path)