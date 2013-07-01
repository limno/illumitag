# Built-in modules #

# Internal modules #

# Third party modules #

# Constants #

###############################################################################
class AnalysisRunner(object):
    """Will run stuff on an aggregate"""

    default_steps = [
        ## Start ###
        {'combine_reads':             {}},
        ## OTUs ###
        {'run_denovo':                {}},
        {'run_open_ref':              {}},
        {'run_progressive':           {}},
        {'run_subsample':             {}},
        ### Plots ###
        #{'make_pool_plots':           {'threads':False}},
        #{'make_outcome_plots':        {'threads':False}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.analysis = parent, parent