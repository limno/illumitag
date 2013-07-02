# Built-in modules #
import os, json
from collections import defaultdict

# Third party modules #
import sh, fastqident
from Bio.SeqIO.FastaIO import FastaWriter

# Internal modules #
from samples import Samples
from outcomes import NoBarcode, OneBarcode, SameBarcode, BadBarcode, GoodBarcode
from illumitag.common import property_cached, AutoPaths
from illumitag.helper.primers import TwoPrimers
from illumitag.fasta.single import FASTA, FASTQ
from illumitag.fasta.paired import PairedFASTQ
from illumitag.fasta.other import QualFile, GroupFile
from illumitag.running.pool_runner import PoolRunner
from illumitag.graphs import pool_plots
from illumitag.helper.barcodes import bar_len

###############################################################################
class Pool(object):
    """An illumina MID is called here a 'pool'."""

    all_paths = """
    /samples/
    /groups/
    /graphs/
    /logs/
    /results/mothur/reads.fasta
    /results/mothur/reads.qual
    /results/mothur/groups.tsv
    /results/qiime/reads.fasta
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.id_name)
    def __str__(self): return self.id_name
    def __iter__(self): return iter(self.children)
    def __len__(self): return self.count
    def __getitem__(self, key): return self.samples[key]

    def __init__(self, json_path, out_dir):
        # Output #
        self.out_dir = out_dir
        # Parse #
        self.json_path = json_path
        with open(json_path) as handle: self.info = json.load(handle)
        # Basic #
        self.account = self.info['uppmax_id']
        self.run_num = self.info['run_num']
        self.run_label = self.info['run_id']
        self.project_short_name = self.info['project']
        self.project_long_name = self.info['project_name']
        # Own attributes #
        self.num = self.info['pool_num']
        self.label = self.info['pool_id']
        self.short_name = self.info['pool']
        self.long_name = self.info['pool_name']
        self.id_name = "run%03d-pool%02d" % (self.run_num, self.num)
        # Special #
        self.samples = Samples(self)
        self.primers = TwoPrimers(self)
        self.loaded = False

    def load(self):
        # Automatic paths #
        self.base_dir = self.out_dir + self.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.samples.load()
        self.primers.load()
        # Raw file pairs #
        self.fwd_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run.label, self.label, self.info['forward_reads'])
        self.rev_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run.label, self.label, self.info['reverse_reads'])
        self.fwd = FASTQ(self.fwd_path)
        self.rev = FASTQ(self.rev_path)
        self.raw = PairedFASTQ(self.fwd.path, self.rev.path, self)
        # Make Outcomes #
        self.no_barcodes   = NoBarcode(self)
        self.one_barcodes  = OneBarcode(self)
        self.same_barcodes = SameBarcode(self)
        self.bad_barcodes  = BadBarcode(self)
        self.good_barcodes = GoodBarcode(self)
        self.outcomes = (self.good_barcodes, self.no_barcodes, self.one_barcodes, self.same_barcodes, self.bad_barcodes)
        self.children = self.outcomes
        # The good reads #
        self.quality_reads = self.good_barcodes.assembled.good_primers.len_filtered
        self.trimmed_reads = self.good_barcodes.assembled.good_primers.trimmed_fastq
        # Primer size #
        self.trim_fwd = bar_len + self.primers.fwd_len
        self.trim_rev = bar_len + self.primers.rev_len
        # Mothur output #
        self.mothur_fasta = FASTA(self.p.mothur_fasta)
        self.mothur_qual = QualFile(self.p.mothur_qual)
        self.mothur_groups = GroupFile(self.p.mothur_groups)
        # Qiime output #
        self.qiime_fasta = FASTA(self.p.qiime_fasta)
        # Runner #
        self.runner = PoolRunner(self)
        # Loaded #
        self.loaded = True

    @property
    def first(self): return self.children[0]

    @property_cached
    def count(self):
        if self.fwd_path.endswith('gz'): return int(sh.zgrep('-c', "^+$", self.fwd_path))
        else: return int(sh.grep('-c', "^+$", self.fwd_path))

    def __call__(self, *args, **kwargs):
        if not self.loaded: self.load()
        self.runner.run(*args, **kwargs)

    def run_slurm(self, *args, **kwargs):
        if not self.loaded: self.load()
        self.runner.run_slurm(*args, **kwargs)

    def create_outcomes(self):
        if not self.loaded: self.load()
        for o in self.outcomes: o.create()
        for r in self.raw.parse_barcodes():
            if len(r.matches) == 0:                              self.no_barcodes.add_pair(r)
            elif len(r.matches) == 1:                            self.one_barcodes.add_pair(r)
            elif r.matches[0].set == r.matches[1].set:           self.same_barcodes.add_pair(r)
            elif r.matches[0].sample is not r.matches[1].sample: self.bad_barcodes.add_pair(r)
            else:                                                self.good_barcodes.add_pair(r)
        for o in self.outcomes: o.close()

    def create_samples(self):
        if not self.loaded: self.load()
        for sample in self.samples: sample.create()
        for r in self.quality_reads.parse_barcodes(): r.first.sample.add_read(r.read)
        for sample in self.samples: sample.close()

    def check_fastq_version(self):
        for o in self.outcomes:
            assert fastqident.detect_encoding(o.fwd_path) == 'sanger'
            assert fastqident.detect_encoding(o.rev_path) == 'sanger'

    def make_mothur_output(self):
        # Trimmed fasta #
        if os.path.exists(self.mothur_fasta.path): os.remove(self.mothur_fasta.path)
        os.symlink(self.trimmed_reads.path, self.mothur_fasta.path)
        # The groups file #
        self.mothur_groups.create()
        for r in self.quality_reads.parse_barcodes():
            sample_name = r.first.sample.short_name
            read_name = '%s\t%s\n' % (r.read.id, sample_name)
            self.mothur_groups.handle.write(read_name)
        self.mothur_groups.close()

    def make_qiime_output(self):
        # Prepare fasta writer #
        handle = open(self.qiime_fasta.path, 'w')
        writer = FastaWriter(handle, wrap=0)
        writer.write_header()
        # Counter #
        counter = defaultdict(int)
        # Do it #
        for r in self.quality_reads.parse_barcodes():
            sample_name = r.first.sample.short_name
            counter[sample_name] += 1
            r.read.id = '%s_%i %s' % (sample_name, counter[sample_name], r.read.id)
            bar_seq = r.read.seq[0:bar_len]
            r.read.description = "orig_bc=%s new_bc=%s bc_diffs=0" % (bar_seq, bar_seq)
            writer.write_record(r.read[self.trim_fwd:-self.trim_rev])
        # Close #
        writer.write_footer()
        handle.close()

    def make_pool_plots(self):
        for cls_name in pool_plots.__all__:
            cls = getattr(pool_plots, cls_name)
            cls(self).plot()