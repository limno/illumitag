# Built-in modules #
import os

# Internal modules #
import illumitag
from illumitag.fasta.single import FASTA

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'
silvamod_path = home + 'glob/16s/silvamod.fasta'
amplified_path = home + 'glob/16s/silvamod_v3_v4.fasta'
aligned_path = home + 'glob/16s/silvamod_v3_v4.align'

# Objects #
silvamod = FASTA(silvamod_path)
amplified = FASTA(amplified_path)
aligned = FASTA(aligned_path)

###############################################################################
def amplify():
    """A function to parse the silvamod 16S database and find the primers within
    the full-length sequences to determine the probable length of our amplified
    region."""
    primers = illumitag.pools[0].primers
    bar_len = illumitag.pools[0].bar_len
    counts = {'success':0, 'only_fwd':0, 'only_rev':0, 'no_primer':0, 'on-the-edge': 0}
    def find_primers(reads):
        for r in reads:
            fwd_match = primers.fwd_regex_uracil.search(str(r.seq))
            rev_match = primers.rev_regex_uracil.search(str(r.seq))
            fwd_pos = fwd_match.start() if fwd_match else None
            rev_pos = rev_match.end() if rev_match else None
            if fwd_pos is not None and rev_pos is not None:
                if fwd_pos<bar_len or rev_pos>len(r)-bar_len:
                    counts['on-the-edge'] += 1
                    continue
                counts['success'] += 1
                yield r[fwd_pos-bar_len:bar_len+rev_pos]
            elif fwd_pos: counts['only_fwd'] += 1
            elif rev_pos: counts['only_rev'] += 1
            else: counts['no_primer'] += 1
    amplified.write(find_primers(silvamod))
    return counts

###############################################################################
def align():
    """A function to align the silvamod 16S database"""
    sh.clustalo('-i', amplified, '-o', aligned)