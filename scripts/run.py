# One project via slurm #
import illumitag; pj = illumitag.projects['test']; pj.run_slurm(time='00:15:00', qos=False)

# Just one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(threads=False)

# Just one function for one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(steps=[{'make_pool_plots':{}}], threads=False)