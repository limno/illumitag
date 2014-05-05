# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
import numpy

# Third party modules #

###############################################################################
class Reporter(object):
    """Reporting statistics on an aggregate."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)
    def __iter__(self): return iter(self.pools)
    def __len__(self): return len(self.pools)
    def __getitem__(self, key): return self.pools[key]

    def __init__(self, aggregate):
        # Attributes #
        self.aggregate, self.parent = aggregate, aggregate
        # Inherited #
        self.pools, self.children = self.parent.pools, self.parent.pools

    @property
    def count(self):
        """The read counts per pool"""
        return sum(map(lambda p: len(p), self.pools))

    @property
    def count_qced(self):
        """The read counts after quality control"""
        return sum(map(lambda p: len(p.quality_reads), self.pools))

    @property
    def avg_quality(self):
        return map(lambda p: p.avg_quality, self.pools)

    @property
    def outcome_percentage(self):
        for o in self.aggregate.outcomes:
            print o.first.doc + ': ' + str(int(round(100*len(o)/self.aggregate.count))) + '%'

    @property
    def avg_assembly_stat(self):
        for bc_outcome in self.aggregate.first.outcomes:
            outcomes = [getattr(p, bc_outcome.short_name) for p in self]
            fail_ratios = [1.0 - len(o.assembled)/len(o) for o in outcomes]
            avg_fail = 100*sum(fail_ratios)/len(fail_ratios)
            print '%s: %.2f %%' % (o.doc, avg_fail)

    @property
    def avg_chimera_test(self):
        for bc_outcome in self.aggregate.first.outcomes:
            outcomes = [getattr(p, bc_outcome.short_name) for p in self]
            uchime_refs =    [o.assembled.good_primers.uchime_ref    for o in outcomes]
            uchime_denovos = [o.assembled.good_primers.uchime_denovo for o in outcomes]
            ref_ratios =    [u.ratio for u in uchime_refs]
            denovo_ratios = [u.ratio for u in uchime_denovos]
            avg_ref    = 100*sum(ref_ratios)/len(ref_ratios)
            avg_denovo = 100*sum(denovo_ratios)/len(denovo_ratios)
            print '%s UCHIME ref: %.2f %%' % (o.doc, avg_ref)
            print '%s UCHIME denovo: %.2f %%' % (o.doc, avg_denovo)

    @property
    def size_fraction_chimeras(self):
        for p in self.pools:
            print "--- Pool %s ---"  % p.short_name
            print "Low refere: ", '%.2g%%' % p.fractions.low.refere.percent
            print "Low denovo: ", '%.2g%%' % p.fractions.low.denovo.percent
            print "Med refere: ", '%.2g%%' % p.fractions.med.refere.percent
            print "Med denovo: ", '%.2g%%' % p.fractions.med.denovo.percent
            print "Big refere: ", '%.2g%%' % p.fractions.big.refere.percent
            print "Big denovo: ", '%.2g%%' % p.fractions.big.denovo.percent

    @property
    def loss_statistics(self):
        for p in self.pools:
            print "--- Pool %s ---"  % p.short_name
            for s in p.loss_statistics: print s.msg % (100 - s.value)
            print 'Assembled over 100: %f%%' % p.good_barcodes.assembled.stats['loss']
        print "--- Average ---"
        for s in ['outcome', 'assembly', 'primers', 'n_filter', 'qual_filter', 'len_filter']:
            print 'Average ' + s + ': %f%%' % (100 - numpy.mean([getattr(p.loss_statistics,s).value for p in self]))
        print 'Assembled over 100: %f%%' % numpy.mean([p.good_barcodes.assembled.stats['loss'] for p in self])
