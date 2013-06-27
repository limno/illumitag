# Built-in modules #
import os, time, getpass, locale

# Third party modules #
import matplotlib

################################################################################
class Graph(object):

    def __init__(self, parent, base_dir=None, short_name=None):
        self.parent = parent
        # Base dir #
        if not base_dir: self.base_dir = self.parent.p.graphs_dir
        else: self.base_dir = base_dir
        # Short name #
        if short_name: self.short_name = short_name
        # Paths #
        self.path = self.base_dir + self.short_name + '.pdf'
        self.csv_path = self.base_dir + self.short_name + '.csv'
        self.json_path = self.base_dir + self.short_name + '.json'

    def save_plot(self, fig, axes, width=18.0, height=10.0, bottom=0.1, top=0.93, left=0.06, right=0.98, sep=()):
        # Adjust #
        fig.set_figwidth(width)
        fig.set_figheight(height)
        fig.subplots_adjust(hspace=0.0, bottom=bottom, top=top, left=left, right=right)
        # Data and source #
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
        # Save it #
        fig.savefig(self.path)