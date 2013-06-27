from illumitag import runs
from illumitag.job import Job
pools = runs[0]
job = Job(pools)
steps = 'a'
job.run()