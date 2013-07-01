# Built-in modules #
import os, sys, time, shutil
from collections import OrderedDict
# Internal modules #
import illumitag
from illumitag.common import get_git_tag
from illumitag.common.tmpstuff import TmpFile

# Third party modules #
import sh

# Constants #

################################################################################
class SLURMCommand(object):
    """Makes launching SLURM commands easy to write and easy to use"""

    script_headers = {
        'bash': "#!/bin/bash -l",
        'python': "#!/usr/bin/env python"
    }

    slurm_header = OrderedDict((
        ('change_dir', {'needed': True,  'tag': '#SBATCH -D %s',          'default': os.path.abspath(os.getcwd())}),
        ('job_name'  , {'needed': False, 'tag': '#SBATCH -J %s',          'default': 'test_slurm'}),
        ('out_file'  , {'needed': True,  'tag': '#SBATCH -o %s',          'default': '/tmp/slurm.out'}),
        ('project'   , {'needed': False, 'tag': '#SBATCH -A %s',          'default': os.environ.get('SLURM_ACCOUNT')}),
        ('time'      , {'needed': True,  'tag': '#SBATCH -t %s',          'default': '0:15:00'}),
        ('machines'  , {'needed': True,  'tag': '#SBATCH -N %s',          'default': '1'}),
        ('cores'     , {'needed': True,  'tag': '#SBATCH -n %s',          'default': '8'}),
        ('email'     , {'needed': False, 'tag': '#SBATCH --mail-user %s', 'default': 'lucas.sinclair@ebc.uu.se'}),
        ('email-when', {'needed': True,  'tag': '#SBATCH --mail-type=%s', 'default': 'END'}),
        ('qos'       , {'needed': False, 'tag': '#SBATCH --qos=%s',       'default': 'short'}),
        ('dependency', {'needed': False, 'tag': '#SBATCH -d %s',          'default': 'after:1'}),
    ))

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, command='print "Hello world"', save_script=False, language='python', **kwargs):
        # Attributes #
        self.save_script = save_script
        self.language = language
        self.name = kwargs.get('job_name')
        # Check command #
        if isinstance(command, list): command = ' '.join(map(str, command))
        # Script header #
        script_header = self.script_headers[language]
        # Slurm header #
        slurm_header = []
        for param, info in self.slurm_header.items():
            if not info['needed'] and not param in kwargs: continue
            if kwargs.get(param): slurm_header.append(info['tag'] % kwargs.get(param))
            else:                 slurm_header.append(info['tag'] % info['default'])
        slurm_header = '\n'.join(slurm_header)
        # Script #
        self.script = '\n'.join((script_header, slurm_header, command))

    def launch(self):
        # Make script #
        if self.save_script:
            path = self.save_script
            with open(self.save_script, 'w') as handle:
                handle.write(self.script)
        else:
            path = TmpFile.from_string(self.script).path
        # Do it #
        sh.sbatch(path, _out=sys.stdout)
        # Clean up #
        if not self.save_script: os.remove(path)

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
        module_dir = os.path.dirname(module.__file__)
        module_name = module.__name__
        repos_dir = os.path.abspath(module_dir + '/../')
        project_name = os.path.basename(repos_dir)
        shutil.copytree(repos_dir, self.log_dir + project_name)
        static_module_dir = self.log_dir + project_name + '/'
        # Archive version #
        self.module_version = module.__version__ + ' ' + get_git_tag(repos_dir)
        # Make script #
        script = """
            import os, sys, time, platform
            print "SLURM: start at {0} on {1}".format(time.asctime(), platform.node())
            sys.path.insert(0, "%s")
            import %s
            print "Using version from {0}".format(os.path.abspath(%s.__file__))
            %s
            print "SBATCH: end at {0}".format(time.asctime())"""
        script = script % (static_module_dir, module_name, module_name, command)
        script = '\n'.join(l.lstrip(' ') for l in script.split('\n') if l)
        script_path = self.log_dir + "run.py"
        # Output #
        output_path = self.log_dir + "run.out"
        # Make SLURM object #
        self.slurm_command = SLURMCommand(script, save_script=script_path, out_file=output_path, **kwargs)

    def launch(self):
        self.slurm_command.launch()