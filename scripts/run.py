# The module #
import illumitag

# One project #
pj = illumitag.projects['test']; pj.run_pools()
# One project via slurm #
pj = illumitag.projects['test']; pj.run_pools_slurm()


# Just one pool #
pj = illumitag.projects['test']; p = pj[0]; p(threads=False)
# Just one pool via slurm #
pj = illumitag.projects['andrea']; p = pj[2]; p.run_slurm()
num = illumitag.projects['inga'].first.run_slurm()
# A few pools #
pj = illumitag.projects['test']; [pool() for pool in pj.pools[1:]]


# Just one function for one pool #
pj = illumitag.projects['test']; p = pj[0]; p(steps=[{'make_pool_plots':{}}], threads=False)
# One function for several pools in parallel #
import playdoh; playdoh.map(lambda p: p.pool_fastqc(), illumitag.projects['evaluation'].pools, cpu=5)
# Just one statistic for one project #
p = illumitag.projects['evaluation']; p.load(); [pl.good_barcodes.relative_std_dev for pl in p]
# Just one graph for one project #
p = illumitag.projects['evaluation']; p.load(); [illumitag.graphs.pool_plots.AssemblyCounts(pl).plot() for pl in p]
# Just one function for one project #
pj = illumitag.projects['evaluation']; pj.load(); [p(steps=[{'check_noalign_counts':{}}]) for p in pj]
# Just one plot for one project #
pj = illumitag.projects['evaluation']; pj.load(); pj.graphs[-1].plot()


# Just one analysis #
pj = illumitag.projects['test']; pj.load(); pj.analysis.run()
pj = illumitag.projects['evaluation']; pj.load(); pj.analysis.run()
# Just one analysis via slurm #
job_id = illumitag.projects['evaluation'].run_analysis_slurm()


# All run graphs #
[r.make_plots() for r in illumitag.runs]


# All pools via slurm #
job_ids = [p.run_slurm() for p in illumitag.pools]
# And analyses via slurm #
ids = [pj.run_analysis_slurm() for pj in illumitag.projects]


# SLURM Report #
illumitag.aggregate.make_slurm_report()


# Regenerate the early exit for one pool #
illumitag.projects['inga'].first.run_slurm([{'make_qiime_output':{}},{'make_mothur_output':{}}])


# Assembly statistics #
ipj = illumitag.projects['evaluation']; pj.load(); pj.reporter.outcome_percentage