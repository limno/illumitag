# Built-in modules #

# Internal modules #
from illumitag.common import flatten

# Third party modules #

################################################################################
names = {
    'Actinobacteria': {
        'species' : ['acI-A1', 'acI-A3', 'acI-A4', 'acI-A5', 'acI-A6', 'acI-A7', 'acI-B1', 'acI-B2', 'acI-B3', 'acI-B4', 'acI-C1', 'acI-C2', 'acTH1-A1', 'acSTL-A1', 'acSTL-A2'],
        'clades' :  ['acSTL', 'acSTL-A', 'acTH1', 'acTH1-A', 'acI', 'acI-A', 'acI-B', 'acI-C']
    },
    'Luna': {
        'species' : ['Luna1-A1', 'Luna1-A2', 'Luna1-A3', 'Luna1-A4', 'Luna2', 'Luna3', 'Myco'],
        'clades' :  ['Luna3', 'Luna3-A', 'acTH2', 'acTH2-A', 'acIII', 'acIII-A', 'Luna1', 'Luna1-A']
    },
    'acV': {
        'species' : ['Iluma-A1', 'Iluma-A2', 'Iluma-B1', 'Iluma-B2', 'Iluma-C1', 'lamia', 'acV-A1', 'acV-A2'],
        'clades' :  ['acV', 'acIV', 'acV-A', 'acIV-D', 'acIV-C', 'acIV-B', 'acIV-A']
    },
    'bacI.': {
        'species' : ['bacI-A1', 'bacI-A2', 'bacI-A3', 'bacI-B1', 'HAL-A1', 'HAL-A2', 'Aquir'],
        'clades' :  ['bacIV-B', 'bacIV-A', 'bacI-B', 'bacI-A', 'bacIV', 'bacI']
    },
    'bacII': {
        'species' : ['Flavo-A1', 'Flavo-A2', 'Flavo-A3', 'Flecto', 'Algor', 'Muci', 'Pedo'],
        'clades' :  ['bacVI', 'bacIII', 'bacV', 'bacII', 'bac-VI-B', 'bac-VI-A', 'bacIII-B', 'bacIII-A', 'bacII-A']
    },
    'alfI': {
        'species' : ['alfI-A1', 'alfI-B1', 'alfI-B2', 'LD12', 'Brev'],
        'clades' :  ['alfII-A', 'alfV-A', 'alfI-B', 'alfI-A', 'alfII', 'alfVIII', 'alfVII', 'alfVI', 'alfV', 'alfI']
    },
    'alfII': {
        'species' : ['Novo-A1', 'Novo-A2', 'Pyxis', 'Sphingo'],
        'clades' :  ['alfIII-A', 'alfIV-B', 'alfIV-A', 'alfIII', 'alfIV']
    },
    'betI': {
        'species' : ['Lhab-A1', 'Lhab-A2', 'Lhab-A3', 'Lhab-A4', 'Rhodo', 'betV-A1'],
        'clades' :  ['betV-A', 'betI-B', 'betl-A', 'betV', 'betI']
    },
    'betII': {
        'species' : ['PnecA', 'PnecB', 'PnecC', 'PnecD', 'LD28', 'betVII-A1', 'Janb', 'betVII-B1', 'betIII-A1', 'betIII-A2'],
        'clades' :  ['betIII-A', 'bet-VII-B', 'betVII-A', 'betIV-A', 'Pnec', 'betVI', 'betIII', 'betVII', 'betIV', 'betII']
    },
    'Gammaproteobacteria': {
        'species' : ['gamII-A1', 'gamII-A2', 'Acin', 'Pseudo-A1', 'Pseudo-A2', 'Steno'],
        'clades' :  ['gamV-A', 'gamIV-A', 'gamIII-A', 'gamII-A', 'gamV', 'gamIV', 'gamIII', 'gamII', 'gamI']
    },
    'Verrucomicrobia': {
        'species' : ['Xip-A1', 'Xip-B1'],
        'clades' :  ['verI-B', 'verI-A', 'verI']
    },
    'Fibrobacteres': {
        'species' : ['CLO-84', 'Fib', 'Cyth', 'Ana'],
        'clades' :  ['Anal-A', 'Cyth-A', 'FibI-A', 'OP10I-A', 'Anal', 'CythI', 'Fib-I', 'OP10I']
    },

}

species_names = flatten([v['species'] for v in names.values()])
clade_names = flatten([v['clades'] for v in names.values()])