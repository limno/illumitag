# Built-in modules #
import sys, time, datetime

# Internal modules #
from illumitag.common import Color

# Third party modules #
import threadpool

# Constants #

###############################################################################
class PoolRunner(object):
    """Will run stuff on a pool"""

    default_steps = [
        ## Outcomes ###
        {'create_outcomes':           {}},
        {'check_fastq_version':       {}},
        ### Assemble ###
        {'assemble':                  {}},
        {'check_noalign_counts':      {}},
        ### Primers ###
        {'flip_reads':                {}},
        {'make_primer_groups':        {}},
        ### Quality ###
        {'discard_reads_with_n':      {}},
        {'quality_filter':            {}},
        {'len_filter':                {}},
        {'trim_barcodes':             {}},
        ### Early exit ##
        {'make_mothur_output':        {}},
        {'make_qiime_output':         {}},
        ### Chimeras ###
        {'check_chimeras':            {}},
        ### FastQC ###
        {'barcode_fastqc':            {}},
        {'assembly_fastqc':           {}},
        ### Plots ###
        {'make_pool_plots':           {'threads':False}},
        {'make_outcome_plots':        {'threads':False}},]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        # Auto color #
        import __main__ as main
        if not hasattr(main, '__file__'): self.color = True
        else: self.color = False

    def run(self, steps=None):
        # Message #
        if self.color: print Color.f_cyn + "Running pool %s" % (self.pool) + Color.end
        else: print "Running pool %s" % self.pool
        # Do steps #
        if not steps: steps = self.default_steps
        for step in steps:
            name, params = step.items()[0]
            fns = self.find_fns(name)
            self.run_step(name, fns, **params)
        # Report success #
        print "Success. Results are in %s" % self.parent.base_dir

    def find_fns(self, name):
        # Functions #
        fns = None
        # Check pool #
        if hasattr(self.pool, name): fns = [getattr(self.pool, name)]
        # Check outcomes #
        elif hasattr(self.pool.first, name): fns = [getattr(o, name) for o in self.pool.outcomes if hasattr(o, name)]
        # Check assemble groups #
        elif hasattr(self.pool.first.first, name): fns = [getattr(ag, name) for o in self.pool.outcomes for ag in o.children if hasattr(ag, name)]
        # Check primer groups #
        elif hasattr(self.pool.first.first.first, name): fns = [getattr(pg, name) for o in self.pool.outcomes for ag in o.children for pg in ag.children if hasattr(pg, name)]
        # Check samples #
        elif hasattr(self.pool.samples.first, name): fns = [getattr(s, name) for s in self.pool.samples]
        # None found #
        if not fns: raise Exception("Could not find function '%s'" % name)
        # Return #
        return fns

    def run_step(self, name, fns, threads=True):
        # Start timer #
        start_time = time.time()
        # Message #
        if self.color: print "Running step: " + Color.f_grn + name + Color.end
        else: print "Running step: " + name
        sys.stdout.flush()
        # Threads #
        if threads:
            self.thpool = threadpool.ThreadPool(8)
            for fn in fns: self.thpool.putRequest(threadpool.WorkRequest(fn))
            self.thpool.wait()
            self.thpool.dismissWorkers(8)
            del self.thpool
        else:
            for fn in fns: fn()
        # Stop timer #
        run_time = datetime.timedelta(seconds=round(time.time() - start_time))
        if self.color: print Color.ylw + "Run time: '%s'" % (run_time) + Color.end
        else: print "Run time: '%s'" % (run_time)
        sys.stdout.flush()
