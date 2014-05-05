# Built-in modules #
import os, time, getpass, locale

# Internal modules #
from illumitag.common.autopaths import FilePath

# Third party modules #
import matplotlib, brewer2mpl

# Constants #
cool_colors = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
cool_colors.reverse()
cool_colors += brewer2mpl.get_map('Set2',    'qualitative', 8).mpl_colors
cool_colors += brewer2mpl.get_map('Set3',    'qualitative', 8).mpl_colors
cool_colors += brewer2mpl.get_map('Pastel1', 'qualitative', 8).mpl_colors
cool_colors += brewer2mpl.get_map('Pastel2', 'qualitative', 8).mpl_colors
cool_colors += brewer2mpl.get_map('Greys',   'sequential', 8).mpl_colors

################################################################################
class Graph(object):
    width = 12.0
    height = 7.0
    bottom = 0.14
    top = 0.93
    left = 0.06
    right = 0.98
    formats = ('pdf')

    def __init__(self, parent, base_dir=None, short_name=None):
        # Save parent #
        self.parent = parent
        # Base dir #
        if not base_dir: self.base_dir = self.parent.p.graphs_dir
        else: self.base_dir = base_dir
        # Short name #
        if short_name: self.short_name = short_name
        # Paths #
        self.path = FilePath(self.base_dir + self.short_name + '.pdf')
        self.csv_path = self.path.replace_extension('csv')
        # Extra #
        self.dev_mode = False

    def save_plot(self, fig, axes, width=None, height=None, bottom=None, top=None, left=None, right=None, sep=()):
        # Attributes or parameters #
        w = width if width != None else self.width
        h = height if height != None else self.height
        b = bottom if bottom != None else self.bottom
        t = top if top != None else self.top
        l = left if left != None else self.left
        r = right if right != None else self.right
        # Adjust #
        fig.set_figwidth(w)
        fig.set_figheight(h)
        fig.subplots_adjust(hspace=0.0, bottom=b, top=t, left=l, right=r)
        # Data and source #
        if self.dev_mode:
            fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
            job_name = os.environ.get('SLURM_JOB_NAME', 'Unnamed')
            user_msg = 'user: %s, job: %s' % (getpass.getuser(), job_name)
            fig.text(0.01, 0.98, user_msg, horizontalalignment='left')
        # Nice digit grouping #
        if 'x' in sep:
            locale.setlocale(locale.LC_ALL, '')
            seperate = lambda x,pos: locale.format("%d", x, grouping=True)
            axes.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
        if 'y' in sep:
            locale.setlocale(locale.LC_ALL, '')
            seperate = lambda x,pos: locale.format("%d", x, grouping=True)
            axes.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(seperate))
        # Save it as different formats #
        for ext in self.formats: fig.savefig(self.path.replace_extension(ext))