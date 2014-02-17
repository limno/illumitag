# Built-in modules #
import os, json, glob, shutil, re

# Internal modules #
import illumitag
from illumitag.common import natural_sort
from illumitag.common.autopaths import AutoPaths, FilePath
from illumitag.common.tmpstuff import TmpFile
from illumitag.fasta.single import FASTQ

# Third party modules #
import sh, pandas, commands
from shell_command import shell_output

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Pyrosample(object):
    """A Pyrosample is a legacy object for the few 454 samples we have"""

    all_paths = """
    /info.json
    /reads.fastq
    """

    def __init__(self, json_path, out_dir):
        # Attributes #
        self.out_dir = out_dir
        self.json_path = FilePath(json_path)
        # Parse #
        with open(json_path) as handle: self.info = json.load(handle)
        # Basic #
        self.run_num = self.info['run_num']
        self.project_short_name = self.info['project']
        self.project_long_name = self.info['project_name']
        # Own attributes #
        self.num = self.info['sample_num']
        self.short_name = self.info['sample']
        self.long_name = self.info['sample_name']
        self.name = 'run%i_sample%i' % (self.run_num, self.num)
        self.group = self.info['group']
        self.id_name = "run%03d-sample%02d" % (self.run_num, self.num)
        # SFF files #
        self.sff_files_info = self.info['files']
        for f in self.sff_files_info:
            if not os.path.exists(f['path']): raise Exception("No file at %s" % f['path'])
        # Automatic paths #
        self.base_dir = self.out_dir + self.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        self.p.info_json.link_from(self.json_path, safe=True)
        # Pool dummy #
        self.pool, self.parent = self, self
        # Other dummy variables #
        self.bar_len = 0
        self.gziped = False
        # FASTQ #
        self.fastq = FASTQ(self.p.reads)

###############################################################################
class Demultiplexer454(object):
    """Will demultiplex a bunch of SFF files"""

    def __init__(self, samples):
        # Save samples #
        self.samples = samples
        # Get all unique SFF files #
        self.sff_paths = set([f['path'] for s in self.samples for f in s.sff_files_info])
        # Make objects of them #
        self.sff_files = [MultiplexedSFF(path) for path in self.sff_paths]
        self.sff_files.sort(key = lambda x: natural_sort(x.name))
        # Two way linked graph #
        for sample in self.samples: sample.sff_links = {}
        for sff in self.sff_files:  sff.sample_links = {}
        # Draw links #
        for sample in self.samples:
            for f in sample.sff_files_info:
                barcode, path = f['mid'], f['path']
                sff_obj = [sff for sff in self.sff_files if sff.path==path][0]
                sample.sff_links[barcode]     = sff_obj
                sff_obj.sample_links[barcode] = sample

    def run(self):
        # Check #
        if glob.glob('454Reads.*.sff'): raise Exception("Demultiplexed files already present in current directory.")
        # Extract #
        for sff in self.sff_files: sff.split()
        # Assign #
        for sample in self.samples: sample.pieces = []
        for sff in self.sff_files: sff.assign_pieces()
        # Report #
        columns = [sff.name for sff in self.sff_files]
        rows = [sample.name for sample in self.samples]
        data = [[sff.count_reads_in(sample) for sff in self.sff_files] for sample in self.samples]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        self.frame['Total'] = self.frame.sum(1)
        # Regroup #
        for sample in self.samples:
            if len(sample.pieces) == 1: shutil.move(sample.pieces[0].path, sample.p.raw_sff)
            else:
                sh.sfffile("-o", sample.p.raw_sff, *[p.path for p in sample.pieces])
                for p in sample.pieces: os.remove(p.path)
        # Save report #
        self.frame.to_csv(self.parent.p.demultiplexing)

###############################################################################
class MultiplexedSFF(FilePath):
    """The MultiplexedSFF object represents an SFF file
    containing several sample pieces distinguishable by
    their barcodes."""

    def __iter__(self): return iter(self.parts)
    def __eq__(self, other): return self.path == other.path

    def __init__(self, path):
        # Basic #
        self.path = path
        self.name = os.path.basename(path)
        # Optional raw output #
        prefix = illumitag.view_dir + 'pyrosamples/raw/' + self.prefix
        self.raw_fasta_path = prefix + ".fasta"
        self.raw_qual_path = prefix + ".qual"
        self.fastq = FASTQ(prefix + ".fastq")

    @property
    def barcode_text(self):
        barcodes = '\n'.join(['mid = "%s", "%s", 0;' % (sample.name,code) for code,sample in self.sample_links.items()])
        return "barcodes_keyword\n{\n%s\n}" % (barcodes)

    def split(self):
        # Call Roche binary #
        barcode_file = TmpFile.from_string(self.barcode_text)
        sh.sfffile("-s", "barcodes_keyword", "-mcf", barcode_file.path, self.path)
        # Check result #
        produced_files = set(sh.glob('454Reads.*.sff'))
        expected_files = set(['454Reads.%s.sff' % (sample.name.upper()) for sample in self.sample_links.values()])
        assert produced_files == expected_files
        # Make piece objects #
        self.pieces = [SamplePiece(p, self) for p in sh.glob('454Reads.*.sff')]
        for piece in self.pieces: piece.rename()
        # Cleanup #
        barcode_file.remove()

    def assign_pieces(self):
        for p in self.pieces: p.sample.pieces.append(p)

    def count_reads_in(self, sample):
        # Find the piece #
        piece = [p for p in self.pieces if p.sample.name == sample.name]
        if not piece: return 0
        else: piece = piece[0]
        # Count reads #
        return int(commands.getstatusoutput('sffinfo -s %s |grep -c ^\>' % piece.path)[1])

    def no_barcode_split(self):
        """Will translate the SFF files without barcodes. Just to explore an SFF file."""
        shell_output('sffinfo -s %s > %s' % (self.path, self.raw_fasta_path))
        shell_output('sffinfo -q %s > %s' % (self.path, self.raw_qual_path))
        sh.fasta_to_fastq(self.raw_fasta_path, self.raw_qual_path, self.fastq)

###############################################################################
class SamplePiece(object):
    """The SamplePiece object represents a piece of a sample extracted
    from a multiplexed file."""

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, path, sff):
        self.path = path
        self.sff = sff
        self.upper_name = re.findall('454Reads\.(.+?)\.sff', self.path)[0]
        self.sample = [s for s in sff.sample_links.values() if s.name.upper()==self.upper_name][0]
        self.name = self.sff.name[:-4] + ":" + self.sample.name

    def rename(self):
        new_path = self.path.replace('454Reads', self.sff.name[:-4])
        new_path = new_path.replace(self.sample.name.upper(), self.sample.name)
        shutil.move(self.path, new_path)
        self.path = new_path