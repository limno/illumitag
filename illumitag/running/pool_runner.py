# Built-in modules #
import os

# Internal modules #
from illumitag.running import Runner
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class PoolRunner(Runner):
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
        {'make_outcome_plots':        {'threads':False}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.pool = parent, parent
        # Default variables #
        self.job_start_time = None
        self.job_end_time = None
        self.job_runtime = None

    def run_slurm(self, steps=None, **kwargs):
        command = """steps = %s
                     pool = [p for p in illumitag.pools if str(p)=='%s'][0]
                     pool(steps)""" % (steps, self)
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = '12:00:00'
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        self.pool.slurm_job = SLURMJob(command, self.pool.p.logs_dir, job_name=str(self.pool), **kwargs)
        return self.pool.slurm_job.launch()

    @property
    def latest_log(self):
        if not self.pool.loaded: self.pool.load()
        def logs():
            for dir_name in os.listdir(self.pool.p.logs_dir):
                dir_path = os.path.join(self.pool.p.logs_dir, dir_name)
                if not os.path.isdir(dir_path): continue
                yield dir_path + '/'
        return max(logs(), key=lambda x: os.stat(x).st_mtime)