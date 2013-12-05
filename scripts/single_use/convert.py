#!/usr/bin/env python2

"""
A script to add some metadata from excel files
to our json files.
"""

# Modules #
import pandas

###############################################################################
#reader = csv.reader(open('metadata_20131123_FINAL.txt'))
#keys = reader.next()[0].split('\t')
#out = [dict(zip(keys, property)) for property in reader]
#with open('result.json', 'w') as handle: json.dumps(handle, out)

# Load data #
df = pandas.io.parsers.read_csv('metadata_domenico.tsv', sep='\t', index_col=0, encoding='utf-8', dtype=str)

# Correspond #
import illumitag
r = illumitag.runs[4]
for pool_num in [1,2,3,4]:
    pool = r[pool_num-1]
    for sample in pool.samples:
        if sample.group_name != 'JDS_2007': sample.extra_metadata = None
        else: sample.extra_metadata = dict(df.loc[sample.short_name])
    with open(pool.id_name, 'w') as handle: handle.write(pool.json)