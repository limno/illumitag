# Built-in modules #

# Internal modules #
from illumitag.clustering.composition import Composition
from illumitag.common.cache import property_cached

# Third party modules #

###############################################################################
class CompositionTips(Composition):
    """The taxa are composed of the lowest level"""

    @property_cached
    def taxa_table(self):
        # Return result #
        return self.frame