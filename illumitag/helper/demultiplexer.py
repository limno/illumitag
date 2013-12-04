"""
Special module to demultiplex the odd hybrid run number 5
"""

# Built-in modules #
import os

# Internal modules #
import illumitag

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Demultiplexer(object):

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent.id_name)

    def __init__(self, fwd, rev, pools):
        self.fwd = FASTQ(fwd)
        self.rev = FASTQ(rev)
        self.pair = FASTQ(rev)
        self.pools = pools

    def run(self):
        pass

###############################################################################
if __name__ == "__main__":
    path = home + "proj35/INBOX/131126_M00485_0087_000000000-A6GWG/Undetermined_indices/Sample_lane1/"
    fwd = "lane1_Undetermined_L001_R1_001.fastq.gz"
    rev = "lane1_Undetermined_L001_R2_001.fastq.gz"
    demulti = Demultiplexer(illumitag.run[5])
    demulti.run()