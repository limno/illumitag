# Built-in modules #
import random
from collections import Counter

# Internal modules #
from illumitag.analysis.otus import OTUs

# Third party modules #
import pandas

# Constants #

#------------------------------------------------------------------------------#
class SubsampledOTUs(OTUs):
    short_name = 'subsampled'
    method = 'Denovo Subsampled OTUs'
    dist_method = 'bray'

    all_paths = OTUs.all_paths + """
    /table/subsampled_float.csv
    """

    def filter_otu_table(self): raise NotImplementedError('')

    def subsample_table(self):
        # Parse #
        otus = pandas.read_csv(self.base_otu.p.table_trimmed, sep = '\t', index_col=0)
        # Drop those below 10 #
        #sums = otus.sum(axis=1)
        #otus = otus.drop(sums[sums < 10].keys())
        # Subsample #
        sums = otus.sum(axis=1)
        down_to = min(sums)
        subotus = pandas.DataFrame(columns=otus.columns, index=otus.index, dtype=int)
        # Do it #
        for sample_name in otus.index:
            row = otus.loc[sample_name]
            weighted_choices = list(row[row != 0].iteritems())
            population = [val for val, count in weighted_choices for i in range(count)]
            sub_pop = random.sample(population, down_to)
            frequencies = Counter(sub_pop)
            new_row = pandas.Series(frequencies.values(), index=frequencies.keys(), dtype=int)
            subotus.loc[sample_name] = new_row
        # Output it #
        subotus.to_csv(self.p.subsampled_table_float, sep='\t', na_rep='0')
        # Cast to integer #
        def lines_as_integer(path):
            handle = open(path)
            yield handle.next()
            for line in handle:
                line = line.split()
                label = line[0]
                values = map(float, line[1:])
                values = map(int, values)
                yield label + '\t' + '\t'.join(map(str,values)) + '\n'
        handle = open(self.p.table_trimmed, 'w')
        handle.writelines(lines_as_integer(self.p.subsampled_table_float))
        handle.close()