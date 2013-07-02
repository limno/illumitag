# One project #
import illumitag; pj = illumitag.projects['test']; pj.run_pools()

# One project via slurm #
import illumitag; pj = illumitag.projects['test']; pj.run_pools_slurm(time='00:15:00')

# A few pools #
import illumitag; pj = illumitag.projects['test']; [pool() for pool in pj.pools[1:]]

# Just one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(threads=False)

# Just one pool via slurm #
import illumitag; pj = illumitag.projects['Andrea']; p = pj[2]; p.run_slurm()
import illumitag; id = illumitag.projects['Inga'].first.run_slurm()

# Just one function for one pool #
import illumitag; pj = illumitag.projects['test']; p = pj[0]; p(steps=[{'make_pool_plots':{}}], threads=False)

# Analysis #
import illumitag; pj = illumitag.projects['test']; pj.load(); pj.analysis.run()

# Analysis via slurm #
import illumitag; pj = illumitag.projects['test']; pj.load(); pj.analysis.run()

# All run type graphs #
import illumitag; [r.make_plots() for r in illumitag.runs]

# All pools via slurm #
import illumitag; job_ids = [p.run_slurm() for p in illumitag.pools]

# All pools and analyses via slurm #
import illumitag
projs = [pj for pj in illumitag.projects if pj.name != 'test']
[pj.load() for pj in projs]
job_ids = sum([pj.run_pools_slurm() for pj in projs], [])
deps = lambda pj: 'afterok:' + ':'.join(str(p.slurm_job.id) for p in pj)
more_ids = [pj.run_analysis_slurm(dependency=deps(pj)) for pj in projs]