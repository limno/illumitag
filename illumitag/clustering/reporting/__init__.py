# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #

# Third party modules #

###############################################################################
class ClusterReporter(object):
    """Reporting statistics on an cluster."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key): return self.samples[key]

    def __init__(self, parent):
        self.cluster, self.parent = parent, parent
        self.samples = self.parent.samples

    @property
    def total_raw_count(self):
        """Number of sequences only after barcode control"""
        return sum([sample.count_raw_reads for sample in self])

    @property
    def total_seq_count(self):
        """Number of sequences after quality control"""
        return sum([len(sample) for sample in self])

    @property
    def otu_count_raw(self):
        """Number of OTUs found before removing undesirables"""
        return len(self.cluster.otu_uparse.centers)

    @property
    def otu_count(self):
        """Number of OTUs found after removing undesirables"""
        return self.cluster.otu_uparse.taxonomy_silva.otu_table.shape

    @property
    def min_max_before_assembly(self):
        """The samples with max and min sequences before assembly"""
        samples = [s for s in self]
        samples.sort(key= lambda x: x.pool.good_barcodes.counter[x.bar_name + 'F'])
        display = lambda x: (x.name, x.short_name, x.pool.good_barcodes.counter[x.bar_name + 'F'])
        return {"Lowest": display(samples[0]),
                "Highest": display(samples[-1])}
