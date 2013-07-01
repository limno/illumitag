# Built-in modules #

# Internal modules #
from illumitag.running import Runner

# Third party modules #

# Constants #

###############################################################################
class AnalysisRunner(Runner):
    """Will run stuff on an aggregate's analysis"""

    default_steps = [
        ### Start ###
        {'combine_reads':             {}},
        ### OTUs ###
        {'run_denovo':                {}},
        #{'run_open_ref':              {}},
        #{'run_progressive':           {}},
        #{'run_subsample':             {}},
        ### Plots ###
        #{'make_pool_plots':           {'threads':False}},
        #{'make_outcome_plots':        {'threads':False}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.analysis = parent, parent