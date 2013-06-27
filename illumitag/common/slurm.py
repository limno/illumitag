# Built-in modules #
import os, copy, time, shutil

# Internal modules #
import illumitag
#from illumitag.common import get_git_tag
from illumitag.common.tmpstuff import TmpFile

# Third party modules #
import pystache, sh

# Constants #

################################################################################
class SLURMCommand(object):
    """Makes launching SLURM commands easy to write and easy to use"""

    script_headers = {
        'bash': "#!/bin/bash -l",
        'python': "#!/usr/bin/env python"
    }

    slurm_header = {
        'change_dir':  {'#SBATCH -D %s':          os.getcwd() + '/'},
        'job_name':    {'#SBATCH -J %s':          'test_slurm'},
        'out_file':    {'#SBATCH -o %s':          '/tmp/slurm.out'},
        'project':     {'#SBATCH -A %s':          os.environ.get('SLURM_ACCOUNT')},
        'time':        {'#SBATCH -t %s':          '0:15:00'},
        'machines':    {'#SBATCH -N %s':          '1'},
        'cores':       {'#SBATCH -n %s':          '1'},
        'email':       {'#SBATCH --mail-user %s': 'lucas.sinclair@ebc.uu.se'},
        'email-when':  {'#SBATCH --mail-type=%s': 'END'},
        'dependency':  {'#SBATCH -d %s':          'after:1'},
    }

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, command='print "Hello world"',
                 save_script=False,
                 language='python',
                 **kwargs):
        # Parameters #
        self.params = copy.copy(self.default_slurm_params)
        for k,v in kwargs.items(): self.params.update[k] = str(v)
        # Script header #
        self.script_header = self.script_headers[self.params['language']]
        # Slurm header #
        slurm_header_template = '\n'.join(l.lstrip(' ') for l in self.header_template.split('\n') if l)
        render = pystache.Renderer(escape=lambda u: u).render
        self.slurm_header = render(slurm_header_template, self.params)
        # Check command #
        if isinstance(self.params['command'], list):
            self.params['command'] = ' '.join(map(str, self.params['command']))
        # Script #
        self.script = '\n'.join((self.script_header, self.slurm_header, self.params['command']))

    def launch(self):
        # Make script #
        if self.params['save_script']:
            path = self.params['save_script']
            with open(path, 'w') as handle:
                handle.write(self.script)
        else:
            path = TmpFile.from_string(self.script).path
        # Do it #
        sh.sbatch(handle.name)
        # Clean up #
        if not self.params['save_script']: os.remove(handle.name)

################################################################################
class SLURMJob(object):
    """Takes care of running a python job through SLURM and logs results.
    Will run it remotely in a new interpreter with a static copy of a module."""

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, command, log_base_dir, module=illumitag, **kwargs):
        # Log directory #
        dir_name = "%4d-%02d-%02d_%02d-%02d-%02d"
        dir_name = dir_name % time.localtime()[0:6]
        self.log_dir = log_base_dir + dir_name + '/'
        os.mkdir(self.log_dir)
        # Copy module there #
        #module_name = module.__name__
        module_dir = os.path.dirname(module.__file__)
        project_dir = os.path.abspath(module_dir + '/../')
        #module_version = module.__version__ + ' ' + get_git_tag(project_dir)
        shutil.copytree(project_dir, self.log_dir)
        # Import that module #

        # Make it #
        self.slurm_command = SLURMCommand(command, **kwargs)

    def launch(self):
        self.slurm_command.launch()