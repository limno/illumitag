# Built-in modules #
import os

# Internal modules #
from illumitag.common import replace_extension
from illumitag.fasta.single import FASTA

# Third party modules #
import sh, shutil

#-----------------------------------------------------------------------------#
class QualFile(FASTA):
    """A single QUAL file somewhere in the filesystem"""
    extension = 'qual'

#-----------------------------------------------------------------------------#
class GroupFile(object):
    """A mothur group file"""

    def __init__(self, path, **kwargs):
        # Save attributes #
        self.path = path
        self.__dict__.update(kwargs)

    def create(self):
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        self.handle = open(self.path, 'w')

    def close(self):
        self.handle.close()

#-----------------------------------------------------------------------------#
class Aligned(FASTA):
    """A single aligned FASTA file somewhere in the filesystem"""

    def __init__(self, aligned_path, orig_path, **kwargs):
        # Save attributes #
        self.aligned_path = aligned_path
        self.path = orig_path
        self.__dict__.update(kwargs)
        # The other paths #
        self.report_path = replace_extension(aligned_path, 'report.txt')
        self.accnos_path = replace_extension(aligned_path, 'accnos.txt')

    def align(self, ref_path):
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=blast, flip=false, processors=8);" % (self.path, ref_path))
        # Move things #
        shutil.move(self.path[:-6] + '.align', self.aligned_path)
        shutil.move(self.path[:-6] + '.align.report', self.report_path)
        shutil.move(self.path[:-6] + '.flip.accnos', self.accnos_path)
        # Clean up #
        if os.path.exists('formatdb.log'): os.remove('formatdb.log')
        if os.path.exists('error.log') and os.path.getsize('error.log') == 0: os.remove('error.log')
        for p in sh.glob('mothur.*.logfile'): os.remove(p)
