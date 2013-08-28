# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
import numpy

# Third party modules #
import playdoh

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
        return sum(map(lambda p: p.count, self.pools))
        return sum(playdoh.map(lambda p: p.count, self.pools, cpu=len(self)))

    @property
    def avg_quality(self):
        return map(lambda p: p.avg_quality, self.pools)
        return playdoh.map(lambda p: p.avg_quality, self.pools, cpu=len(self))

    @property
    def outcome_percentage(self):
        for o in self.aggregate.outcomes:
            print o.first.doc + ': ' + str(int(round(100*len(o)/self.aggregate.count))) + '%'
        #assert sum(map(len, p.outcomes)) == p.count

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
