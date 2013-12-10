#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to convert some custom excel files to our json file format.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./xls_to_json.py < metadata.tsv > metadata.json
"""

# Modules #
import pandas, simplejson, sys, codecs

###############################################################################
#df = pandas.io.parsers.read_csv(sys.stdin, sep='\t', index_col=0, encoding='mac_roman')
#result = {i: dict(row) for i,row in df.iterrows()}
#sys.stdout.write(json.dumps(result, indent=4, ensure_ascii=False))

###############################################################################
units = {
        "Temperature": ['temperature', 'Celsius'],
        "Oxygen (O2)": ['oxygen', 'mg/L', 'O2'],
        "Chla": ['Chla', 'µg/l'],
        "Na": ['Sodium', 'mg/l', 'Na'],
        "Wind direction": ['wind direction', 'cardinal'],
        "TOC": ['TOC', 'mg/l'],
        "Ptot": ['Ptot', 'µg/l'],
        "pH": ['pH', '-log10([H+])'],
        "TSS": ['TSS', 'g/L', '±50'],
        "Oxygen (O2).1": ['oxygen_sat', '%', '±0.2'],
        "HCO3": ['temperature', 'mg/l', 'HCO3'],
        "Wind strength": ['wind strength', 'beaufort scale'],
        "Hardness": ['hardness', '°dH'],
        "Mg": ['magnesium', 'mg/l', 'Mg'],
        "Clouds": ['clouds', '% coverage'],
        "K": ['potassium', 'mg/l', 'K'],
        "Secchi depth": ['secchi depth', 'cm'],
        "SO4": ['sulfate', 'mg/l', 'SO4'],
        "CO3": ['carbonate', 'mg/l', 'CO3'],
        "Cl": ['chlorine', 'mg/l', 'Cl'],
        "Color": ['color', 'mg/l'],
        "Ca": ['calcium', 'mg/l', 'Ca'],
        "Num Measure": ['name', 'de'],
        "Unnamed: 12": ['GPS E', 'E'],
        "Conductivity": ['conductance', 'μS/cm'],
        "GPS Location": ['GPS N', 'N']}

###############################################################################
# Parse #
df = pandas.io.parsers.read_csv('metadata.tsv', sep='\t', index_col=0, encoding='mac_roman')
data = {i: dict(row) for i,row in df.iterrows()}

# Format #
result = {}
for k,row in data.items():
    result[k] = {}
    for measure, value in row.items():
        info = units[measure]
        result[k][info[0]] = [value, info[1]]

# Write #
with open('metadata.json' ,'w') as handle:
    handle = codecs.getwriter('utf8')(handle)
    simplejson.dump(result, handle, indent=4, ensure_ascii=False, encoding="utf-8")