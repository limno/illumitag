# Built-in modules #
import sys, time, datetime

# Internal modules #
from common import Color

# Third party modules #
import threadpool

# Constants #

###############################################################################
class PoolRunner(object):
    """Will run stuff on a pool"""

    default_steps = [
        ### Barcodes ###
        {'create_groups':             {}},
        {'check_fastq_version':       {}},
        {'barcode_stats':             {}},
        {'plot_barcode_barstack':     {'threads':False}},
        {'plot_barcode_hist':         {'threads':False}},
        {'plot_salvage_hist':         {'threads':False}},
        {'plot_missmatch_reg':        {'threads':False}},
        {'barcode_fastqc':            {}},
        ### Assemble ###
        {'assemble':                  {}},
        {'check_noalign_counts':      {}},
        {'plot_assembly_counts':      {'threads':False}},
        {'plot_assembly_distrib':     {'threads':False}},
        {'assembly_fastqc':           {}},
        ### Primers ###
        {'flip_reads':                {}},
        {'make_primer_groups':        {}},
        {'plot_ass_primer_pos':       {'threads':False}},
        {'plot_unass_primer_pos':     {'threads':False}},
        {'plot_primer_counts':        {'threads':False}},
        ### Quality ###
        {'discard_reads_with_n':      {}},
        {'quality_filter':            {}},
        {'len_filter':                {}},
        {'plot_reads_with_n':         {'threads':False}},
        {'plot_quality_filter':       {'threads':False}},
        {'plot_len_filter':           {'threads':False}},
        ### Early exit ##
        {'make_mothur_output':        {}},
        {'make_qiime_output':         {}},
        {'plot_reads_that_pass_hist': {'threads':False}},
        ### Chimeras ###
        {'check_chimeras_ref':        {'threads':False}},
        {'plot_chimeras_ref':         {'threads':False}},
        {'check_chimeras_denovo':     {'threads':False}},
        {'plot_chimeras_denovo':      {'threads':False}},]

    def __init__(self, pool):
        self.pool = pool

    def run(self, steps=None):
        if not steps: steps = self.default_steps
        for step in steps:
            name, params = step.items()[0]
            fns = self.find_fns(name)
            self.run_step(name, fns, **params)

    def find_fns(self, name):
        # Functions #
        fns = None
        # Check pool #
        if hasattr(self.pool, name): fns = [getattr(self.pool, name)]
        # Check outcomes #
        elif hasattr(self.first, name): fns = [getattr(o, name) for o in self.pool.outcomes if hasattr(o, name)]
        # Check assemble groups #
        elif hasattr(self.first.first, name): fns = [getattr(ag, name) for o in self.pool.outcomes for ag in o.children if hasattr(ag, name)]
        # Check primer groups #
        elif hasattr(self.first.first.first, name): fns = [getattr(pg, name) for o in self.pool.outcomes for ag in o.children for pg in ag.children if hasattr(pg, name)]
        # None found #
        if not fns: raise Exception("Could not find function '%s'" % name)

    def run_step(self, name, fns, threads=False, color=False):
        # Start timer #
        start_time = time.time()
        # Auto color #
        import __main__ as main
        if not hasattr(main, '__file__'): color = True
        # Message #
        if color: print "Running step: " + Color.f_grn + name + Color.end
        else: print "Running step: " + name
        sys.stdout.flush()
        # Threads #
        if threads:
            self.pool = threadpool.ThreadPool(8)
            for fn in fns: self.pool.putRequest(threadpool.WorkRequest(fn))
            self.pool.wait()
            self.pool.dismissWorkers(8)
            del self.pool
        else:
            for fn in fns: fn()
        # Stop timer #
        run_time = datetime.timedelta(seconds=round(time.time() - start_time))
        if color: print Color.ylw + "Run time: '%s'" % (run_time) + Color.end
        else: print "Run time: '%s'" % (run_time)
        sys.stdout.flush()
