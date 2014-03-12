# Built-in modules #

# Internal modules #
import illumitag

# Third party modules #

###############################################################################
# Pyro cluster #
samples = [s for s in illumitag.pyrosamples]
samples += illumitag.runs[1][0][0:8]
samples += illumitag.runs[1][1][0:8]
samples += illumitag.runs[2][0][0:8]
samples += illumitag.runs[3][7][30:33]
samples += illumitag.runs[3][7][13:16]
samples += illumitag.runs[3][7][43:46]
pyro_comparison = illumitag.clustering.Cluster(samples, 'pyro_comparison')

# Domenico #
samples = [s for pool in illumitag.runs[4][0:3] for s in pool.samples if s.used]
samples += illumitag.runs[4][3][0:13]
domenico = illumitag.clustering.Cluster(samples, 'domenico')

# Mixed eval #
#samples = [s for s in illumitag.presamples]
#samples += illumitag.runs[1][0][0:8]
#samples += illumitag.runs[1][1][0:8]
#samples += illumitag.runs[2][0][0:8]
#mixed_evaluation = illumitag.clustering.Cluster(samples, 'mixed_evaluation')

# Other clusters #
#new_lab_test_with = illumitag.clustering.Cluster(illumitag.presamples, 'new_lab_test_with')

# Soda lakes #
temporal = [s for s in illumitag.runs[3][7].samples if s.used]
spatial = [s for s in illumitag.runs[4][3].samples if s.used and s.group_name != 'JDS_2007']
samples = temporal + spatial
soda = illumitag.clustering.Cluster(samples, 'soda')

# Inga's cluster #
samples = [s for s in illumitag.runs[3][6].samples if s.used]
inga = illumitag.clustering.Cluster(samples, 'inga')

# Anna's cluster #
samples = [s for s in illumitag.runs[4][4].samples if s.used]
samples += [s for s in illumitag.runs[4][5].samples if s.used]
anna = illumitag.clustering.Cluster(samples, 'anna')

# Jerome's cluster #
samples = [s for s in illumitag.runs[4][6].samples if s.used]
samples += [s for s in illumitag.runs[4][7].samples if s.used]
jerome = illumitag.clustering.Cluster(samples, 'jerome')

# Monica's cluster #
samples =  [s for s in illumitag.runs[5][3].samples if s.used]
samples += [s for s in illumitag.runs[5][4].samples if s.used]
samples += [s for s in illumitag.runs[5][5][0:11] if s.used]
monica = illumitag.clustering.Cluster(samples, 'monica')
