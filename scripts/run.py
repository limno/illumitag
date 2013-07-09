# One project #
import illumitag; pj = illumitag.projects['test']; pj.run_pools()
# One project via slurm #
import illumitag; pj = illumitag.projects['test']; pj.run_pools_slurm()


# Just one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(threads=False)
# Just one pool via slurm #
import illumitag; pj = illumitag.projects['andrea']; p = pj[2]; p.run_slurm()
import illumitag; num = illumitag.projects['inga'].first.run_slurm()
# A few pools #
import illumitag; pj = illumitag.projects['test']; [pool() for pool in pj.pools[1:]]


# Just one function for one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(steps=[{'make_pool_plots':{}}], threads=False)


# Just one analysis #
import illumitag; pj = illumitag.projects['test']; pj.load(); pj.analysis.run()
import illumitag; pj = illumitag.projects['evaluation']; pj.load(); pj.analysis.run()
# Just one analysis via slurm #
import illumitag; job_id = illumitag.projects['evaluation'].run_analysis_slurm()


# All run graphs #
import illumitag; [r.make_plots() for r in illumitag.runs]


# All pools via slurm #
import illumitag; job_ids = [p.run_slurm() for p in illumitag.pools]
# And analyses via slurm #
ids = [pj.run_analysis_slurm() for pj in illumitag.projects]


# SLURM Report #
import illumitag; illumitag.aggregate.make_slurm_report()