# Built-in modules #
import os, json

# Internal modules #
from illumitag.groups.samples import Samples
from illumitag.fasta.paired import PairedFASTQ
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.fasta.single import FASTQ

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Pyrosample(object):
    """A Pyrosample is a legacy object for the few 454 samples we have"""

    all_paths = """
    /info.json
    /reads.fastq
    """

    def __init__(self, json_path, out_dir):
        # Attributes #
        self.out_dir = out_dir
        self.json_path = FilePath(json_path)
        # Parse #
        with open(json_path) as handle: self.info = json.load(handle)
        # Basic #
        self.account = self.info['uppmax_id']
        self.run_num = self.info['run_num']
        self.run_label = self.info['run_id']
        self.project_short_name = self.info['project']
        self.project_long_name = self.info['project_name']
        self.file_name = self.info['reads']
        # Own attributes #
        self.num = self.info['sample_num']
        self.label = self.info['sample_id']
        self.short_name = self.info['sample']
        self.long_name = self.info['sample_name']
        self.name = 'run%i_sample%i' % (self.run_num, self.num)
        self.group = self.info['group']
        self.id_name = "run%03d-sample%02d" % (self.run_num, self.num)
        # Automatic paths #
        self.base_dir = self.out_dir + self.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        self.p.info_json.link_from(self.json_path, safe=True)
        # Samples dummy #
        self.info['samples'] = [{"name":self.short_name, "used":1, "group":self.group,
                                 "dummy":1, "num":self.num, "fwd":"", "rev":""}]
        self.samples = Samples(self)
        self.samples.load()
        # Pool dummy #
        self.pool, self.parent = self, self
        # Barcode length #
        self.bar_len = 0
        # Files #
        self.sff_path = home + "ILLUMITAG/INBOX/%s/%s/%s" % (self.run_label, self.label, self.fwd_name)
        self.rev_path = home + "ILLUMITAG/INBOX/%s/%s/%s" % (self.run_label, self.label, self.rev_name)
        self.gziped = False
        self.sff = SFF(self.p.fwd)
        self.rev = FASTQ(self.p.rev)
        # FASTQ #
        self.fastq = PairedFASTQ(self.fwd.path, self.rev.path, self)

###############################################################################
class SFF(object):
    """A SFF file somewhere on the file system"""

    all_paths = """
    /info.json
    /reads.fastq
    """

    def __init__(self, json_path, out_dir):
        # Attributes #
        self.out_dir = out_dir
