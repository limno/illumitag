# Built-in modules #
import json
from collections import defaultdict

# Third party modules #
import sh, fastqident
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

# Internal modules #
from common import property_cached, AutoPaths
from groups.samples import Samples
from group.outcomes import NoBarcode, OneBarcode, SameBarcode, BadBarcode, GoodBarcode
from helper.primers import TwoPrimers
from fasta.single import FASTA, FASTQ
from fasta.paired import PairedFASTQ
from fasta.other import QualFile, GroupFile
from runner import PoolRunner
from graphs import pool_plots

###############################################################################
class Pool(object):
    """An illumina MID is called here a 'pool'."""

    all_paths = """
    /samples/
    /groups/
    /graphs/
    /results/mothur/reads.fasta
    /results/mothur/reads.qual
    /results/mothur/groups.tsv
    /results/qiime/reads.fasta
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.id_name)
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
        self.id_name = "run%i-pool%i" % (self.num, self.run_num)
        # Special #
        self.samples = Samples(self)
        self.primers = TwoPrimers(self)

    def load(self):
        # Children #
        self.samples.load()
        self.primers.load()
        # Raw file pairs #
        self.fwd_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run.label, self.label, self.info['forward_reads'])
        self.rev_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run.label, self.label, self.info['reverse_reads'])
        self.fwd = FASTQ(self.fwd_path)
        self.rev = FASTQ(self.rev_path)
        self.raw = PairedFASTQ(self.fwd.path, self.rev.path, self)
        # Automatic paths #
        self.base_dir = self.out_dir + self.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.fwd = FASTQ(self.fwd_path)
        self.rev = FASTQ(self.rev_path)
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
        # Mothur output #
        self.mothur_fasta = FASTA(self.p.mothur_fasta)
        self.mothur_qual = QualFile(self.p.mothur_qual)
        self.mothur_groups = GroupFile(self.p.mothur_groups)
        # Qiime output #
        self.qiime_fasta = FASTA(self.p.qiime_fasta)
        # Runner #
        self.runner = PoolRunner(self)

    @property
    def first(self): return self.children[0]

    @property_cached
    def count(self):
        if self.fwd_path.endswith('gz'): return int(sh.zgrep('-c', "^+$", self.fwd_path))
        else: return int(sh.grep('-c', "^+$", self.fwd_path))

    def create_outcomes(self):
        for o in self.outcomes: o.create()
        for r in self.raw.parse_barcodes():
            if len(r.matches) == 0:                        self.no_barcodes.add_pair(r)
            elif len(r.matches) == 1:                      self.one_barcodes.add_pair(r)
            elif r.matches[0].set == r.matches[1].set:     self.same_barcodes.add_pair(r)
            elif r.matches[0].index != r.matches[1].index: self.bad_barcodes.add_pair(r)
            else:                                          self.good_barcodes.add_pair(r)
        for o in self.outcomes: o.close()

    def create_samples(self):
        for sample in self.samples: sample.create()
        for r in self.quality_reads.parse_barcodes():
            self.samples[r.first.index].add_read(r)
        for sample in self.samples: sample.close()

    def check_fastq_version(self):
        for group in self.groups:
            assert fastqident.detect_encoding(group.fwd_path) == 'sanger'
            assert fastqident.detect_encoding(group.rev_path) == 'sanger'

    def make_mothur_output(self):
        # Trimmed fasta #
        self.mothur_fasta.create()
        for r in self.quality_reads:
            SeqIO.write(r[self.trim_start:-self.trim_end], self.mothur_fasta.handle, 'fasta')
        self.mothur_fasta.close()
        # Trimmed qual #
        self.mothur_qual.create()
        for r in self.quality_reads:
            SeqIO.write(r[self.trim_start:-self.trim_end], self.mothur_qual.handle, 'qual')
        self.mothur_qual.close()
        # The groups file #
        self.mothur_groups.create()
        for r in self.quality_reads.parse_barcodes():
            bar_name = self.barcodes.names[r.first.index]
            read_name = '%s\tpool%i-%s\n' % (r.read.id, self.num, bar_name)
            #read_name = '%s\t%s-%s\n' % (r.read.id, self.label_name, self.sample_labels[bar_name])
            self.mothur_groups.handle.write(read_name)
        self.mothur_groups.close()

    def make_qiime_output(self):
        # Prepare fasta writer #
        fasta_handle = open(self.qiime_fasta.path, 'w')
        fasta = FastaWriter(fasta_handle, wrap=0)
        fasta.write_header()
        # Counter #
        counter = defaultdict(int)
        # Do it #
        for r in self.quality_reads.parse_barcodes():
            bar_name = self.barcodes.names[r.first.index]
            counter[bar_name] += 1
            r.read.id = 'pool%i-%s_%i %s' % (self.num, bar_name, counter[bar_name], r.read.id)
            #r.read.id = '%s-%s_%i %s' % (self.label_name, self.sample_labels[bar_name], counter[bar_name], r.read.id)
            r.read.description = "orig_bc=%s new_bc=%s bc_diffs=0" % (r.first.bar, r.first.bar)
            fasta.write_record(r.read[self.trim_start:-self.trim_end])
        # Close #
        fasta.write_footer()
        fasta_handle.close()

    def make_pool_plots(self):
        for cls_name in pool_plots.__all__:
            cls = getattr(pool_plots, cls_name)
            cls(self).plot()