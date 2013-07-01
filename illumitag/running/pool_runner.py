# Built-in modules #

# Internal modules #
from illumitag.running import Runner

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
