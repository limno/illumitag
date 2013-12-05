# Built-in modules #

# Internal modules #
from illumitag.running import Runner
from illumitag.common.slurm import SLURMJob

# Third party modules #

# Constants #

###############################################################################
class PresampleRunner(Runner):
    """Will run stuff on a Presample"""
    default_time = '4:00:00'

    default_steps = [
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
        ### Chimeras ###
        {'check_chimeras':            {}},
        ### FastQC ###
        {'barcode_fastqc':            {}},
        {'assembly_fastqc':           {}},
        ### Plots ###
        {'make_presample_plots':      {'threads':False}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.presample = parent, parent
        # Default variables #
        self.job_start_time = None
        self.job_end_time = None
        self.job_runtime = None

    def find_fns(self, name):
        # Check quality reads #
        if hasattr(self.parent.quality_reads, name): return [getattr(self.parent.quality_reads, name)]
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
                     presample = [p for p in illumitag.presample if str(p)=='%s'][0]
                     presample(steps)""" % (steps, self.parent)
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'dependency' not in kwargs: kwargs['dependency'] = 'singleton'
        job_name = "illumitag_%s" % self.parent
        self.parent.slurm_job = SLURMJob(command, self.parent.p.logs_dir, job_name=job_name, **kwargs)
        return self.parent.slurm_job.run()

    def empty_logs(self):
        self.parent.p.logs_dir.remove()
        self.parent.p.logs_dir.create()