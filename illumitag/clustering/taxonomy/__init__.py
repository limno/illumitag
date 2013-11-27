# Built-in modules #

# Internal modules #
from illumitag.common.autopaths import AutoPaths
from illumitag.fasta.single import FASTA

# Third party modules #

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
        return dict((k,v[level]) for k,v in self.assignments.items() if len(v) > level)
