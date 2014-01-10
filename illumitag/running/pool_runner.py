# Built-in modules #

# Internal modules #
from illumitag.running import Runner
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class PoolRunner(Runner):
    """Will run stuff on a pool"""
    default_time = '12:00:00'

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
        {'length_filter':             {}},
        {'trim_barcodes':             {}},
        ### Early exit ##
        {'filter_unused':             {}},
        {'trim_primers':              {}},
        {'make_mothur_output':        {}},
        {'make_qiime_output':         {}},
        ### Samples ###
        {'create_samples':            {}},
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

    def find_fns(self, name):
        # Check quality reads #
        if hasattr(self.parent.quality_reads, name): return [getattr(self.parent.quality_reads, name)]
        # Check samples #
        elif hasattr(self.parent.samples.first, name): return [getattr(s, name) for s in self.pool.samples]
        # Super #
        return Runner.find_fns(self, name)

    def run_slurm(self, steps=None, **kwargs):
        # Test case #
        if self.parent.project.name == 'test':
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Script #
        command = """steps = %s
                     pool = [p for p in illumitag.pools if str(p)=='%s'][0]
                     pool(steps)""" % (steps, self.parent)
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        job_name = "illumitag_%s" % self.pool
        self.parent.slurm_job = SLURMJob(command, self.parent.p.logs_dir, job_name=job_name, **kwargs)
        return self.parent.slurm_job.run()

    def empty_logs(self):
        self.parent.p.logs_dir.remove()
        self.parent.p.logs_dir.create()