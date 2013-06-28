# Built-in modules #
import re, os, time, shutil

# Internal modules #
from illumitag.common import find_files_by_regex, prepend_to_file, append_to_file

# Third party modules #

# Constants #
primer_silva_start = 6427
primer_silva_end = 23440

# Functions #
def process_log_file(name, directory, start_time=None):
    # Find one log file #
    if not directory.endswith('/') : directory = directory + '/'
    logs = list(find_files_by_regex('mothur.[0-9]+.logfile'))
    if len(logs) > 1: raise Exception("Found more than one mothur log file in '%s'" % os.getcwd())
    path = logs[0]
    # Add time #
    prepend_to_file(path, "Start time: " + start_time + "\n")
    append_to_file(path, "\nEnd time: " + time.asctime())
    # Move it #
    destination = directory + name
    shutil.move(path, directory + name)
    # Check for error #
    raw = open(destination, 'r').read()
    if re.findall('\[ERROR\]', raw, re.M):
        shutil.move(destination, directory + name.upper())
        raise Exception("Mothur log file reports an error in step '%s'" % name)
    if os.path.exists(directory + name.upper()): os.remove(directory + name.upper())