# Built-in modules #
import os, sys, time, datetime

# Internal modules #
from illumitag.common.color import Color

# Third party modules #
import threadpool

# Constants #

###############################################################################
class Runner(object):
    """General purpose runner"""

    @property
    def color(self):
        import __main__ as main
        if not hasattr(main, '__file__'): return True
        return False

    def run(self, steps=None, **kwargs):
        # Message #
        if self.color: print Color.f_cyn + "Running %s" % (self.parent) + Color.end
        else: print "Running %s" % self.parent
        # Do steps #
        if not steps: steps = self.default_steps
        for step in steps:
            name, params = step.items()[0]
            params.update(kwargs)
            fns = self.find_fns(name)
            self.run_step(name, fns, **params)
        # Report success #
        print "Success. Results are in %s" % self.parent.base_dir

    def find_fns(self, name):
        # Functions #
        fns = None
        # Check pool #
        if hasattr(self.parent, name): fns = [getattr(self.parent, name)]
        # Check outcomes #
        elif hasattr(self.parent.first, name): fns = [getattr(o, name) for o in self.parent.children if hasattr(o, name)]
        # Check assemble groups #
        elif hasattr(self.parent.first.first, name): fns = [getattr(ag, name) for o in self.parent.children for ag in o.children if hasattr(ag, name)]
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
        if threads and len(fns) > 1:
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

    @property
    def latest_log(self):
        if not self.parent.loaded: self.parent.load()
        def logs():
            for dir_name in os.listdir(self.parent.p.logs_dir):
                dir_path = os.path.join(self.parent.p.logs_dir, dir_name)
                if not os.path.isdir(dir_path): continue
                yield dir_path + '/'
        return max(logs(), key=lambda x: os.stat(x).st_mtime)