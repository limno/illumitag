# Futures #
from __future__ import division

# Built-in modules #
import os, re, glob, random, collections

# Third party modules #
import sh, numpy

# Expositions #
from color import Color
from autopaths import AutoPaths
from tmpstuff import TmpFile

################################################################################
def get_git_tag(directory):
    if os.path.exists(directory + '/.git'):
        return sh.git("--git-dir=" + directory + '/.git', "describe", "--tags", "--dirty", "--always").strip('\n')
    else:
        return None

################################################################################
def reverse_compl_with_name(old_seq):
    """Reverse a sequence, but keep its name intact"""
    new_seq = old_seq.reverse_complement()
    new_seq.id = old_seq.id
    new_seq.description = old_seq.description
    return new_seq

################################################################################
def downsample_freq(freq, k):
    """Downsample a collections of things through their frequency"""
    pop = flatten([[key]*val for key,val in freq.items()])
    smaller_pop = random.sample(pop, k)
    return collections.Counter(smaller_pop)

################################################################################
def isubsample(full_sample, k, full_sample_len=None):
    """Downsample an enumerable list of things"""
    # Determine length #
    if not full_sample_len: full_sample_len = len(full_sample)
    # Check size coherence #
    if not 0 <= k <= full_sample_len:
        raise ValueError('Required that 0 <= k <= full_sample_length')
    # Do it #
    picked = 0
    for i, element in enumerate(full_sample):
        prob = (k-picked) / (full_sample_len-i)
        if random.random() < prob:
            yield element
            picked += 1
    # Did we pick the right amount #
    assert picked == k

################################################################################
def imean(numbers):
    """Iterative mean"""
    count = 0
    total = 0
    for num in numbers:
        count += 1
        total += num
    return float(total)/count

################################################################################
def find_files_by_regex(regex, directory='./'):
    """Search one directory for all files matching a pattern"""
    for name in glob.glob(directory + "*"):
        if re.match(regex, os.path.basename(name)): yield name

###############################################################################
def replace_extension(path, new_ext):
    if not new_ext.startswith('.'): new_ext = '.' + new_ext
    base, ext = os.path.splitext(path)
    return base + new_ext

################################################################################
flatten = lambda x: [item for sublist in x for item in sublist]

################################################################################
def moving_average(interval, windowsize):
    window = numpy.ones(int(windowsize))/float(windowsize)
    return numpy.convolve(interval, window, 'valid')

################################################################################
def prepend_to_file(path, data, bufsize=1<<15):
    # Backup the file #
    backupname = path + os.extsep + 'bak'
    # Remove previous backup if it exists #
    try: os.unlink(backupname)
    except OSError: pass
    os.rename(path, backupname)
    # Open input/output files,  note: outputfile's permissions lost #
    with open(backupname) as inputfile:
        with open(path, 'w') as outputfile:
            outputfile.write(data)
            buf = inputfile.read(bufsize)
            while buf:
                outputfile.write(buf)
                buf = inputfile.read(bufsize)
    # Remove backup on success #
    os.remove(backupname)

def append_to_file(path, data):
    with open(path, "a") as handle:
        handle.write(data)

################################################################################
def tail(path, window=20):
    with open(path, 'r') as f:
        BUFSIZ = 1024
        f.seek(0, 2)
        num_bytes = f.tell()
        size = window + 1
        block = -1
        data = []
        while size > 0 and num_bytes > 0:
            if num_bytes - BUFSIZ > 0:
                # Seek back one whole BUFSIZ
                f.seek(block * BUFSIZ, 2)
                # Read BUFFER
                data.insert(0, f.read(BUFSIZ))
            else:
                # File too small, start from beginning
                f.seek(0,0)
                # Only read what was not read
                data.insert(0, f.read(num_bytes))
            linesFound = data[0].count('\n')
            size -= linesFound
            num_bytes -= BUFSIZ
            block -= 1
        return '\n'.join(''.join(data).splitlines()[-window:])

def head(path, window=20):
    with open(path, 'r') as handle:
        return ''.join(handle.next() for line in xrange(window))

################################################################################
def property_cached(f):
    """Decorator for properties evaluated only once.
    It can be used to created a cached property like this::

        class Employee(object):
            @property_cached
            def salary(self):
                return 8000

        bob = Employee()
        print bob.salary
    """
    def get_method(self):
        try:
            return self.__cache__[f.func_name]
        except AttributeError:
            self.__cache__ = {}
            x = self.__cache__[f.func_name] = f(self)
            return x
        except KeyError:
            x = self.__cache__[f.func_name] = f(self)
            return x
    get_method.__doc__ = f.__doc__
    return property(get_method)

###############################################################################
def split_thousands(s, tSep='\'', dSep='.'):
    """
    Splits a number on thousands.

    >>> split_thousands(1000012)
    "1'000'012"
    """
    if s is None: return 0
    if not isinstance(s, str): s = str(s)
    cnt=0
    numChars=dSep+'0123456789'
    ls = len(s)
    while cnt < ls and s[cnt] not in numChars: cnt += 1
    lhs = s[0:cnt]
    s = s[cnt:]
    if dSep == '':
        cnt = -1
    else:
        cnt = s.rfind(dSep)
    if cnt > 0:
        rhs = dSep + s[cnt+1:]
        s = s[:cnt]
    else:
        rhs = ''
    splt=''
    while s != '':
        splt= s[-3:] + tSep + splt
        s = s[:-3]
    return lhs + splt[:-1] + rhs
