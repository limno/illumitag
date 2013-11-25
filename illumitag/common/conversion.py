# Built-in modules #

# Internal modules #

# Third party modules #
import pandas, numpy

################################################################################
def r_matrix_to_dataframe(matrix):
    cols = list(matrix.colnames)
    rows = list(matrix.rownames)
    return pandas.DataFrame(numpy.array(matrix), index=rows, columns=cols)