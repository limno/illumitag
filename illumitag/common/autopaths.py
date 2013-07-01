# Built-in modules #
import os, tempfile, re

################################################################################
class AutoPaths(object):
    """
    You can use this class like this when making pipelines:

        class Sample(object):
            all_paths = '''
                /raw/raw.sff
                /raw/raw.fastq
                /clean/trim.fastq
                /clean/clean.fastq'''

            def __init__(self, base_dir):
                self.p = AutoPaths(base_dir, self.all_paths)

            def clean(self):
                shutil.move(self.p.raw_sff, self.p.clean_fastq)
    """

    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self._base_dir)

    def __init__(self, base_dir, all_paths):
        # Attributes #
        self._base_dir = base_dir
        self._all_paths = all_paths
        self._tmp_dir = tempfile.gettempdir() + '/'
        # Parse input #
        self._paths = [PathItems(p.lstrip(' '),base_dir) for p in all_paths.split('\n')]

    def __getattr__(self, key):
        # Let magic pass through #
        if key.startswith('__') and key.endswith('__'): return object.__getattr__(key)
        # Special cases #
        if key.startswith('_'): return self.__dict__[key]
        # Temporary items #
        if key == 'tmp_dir': return self.__dict__['_tmp_dir']
        if key == 'tmp': return self.__dict__['tmp']
        # Search #
        items = key.split('_')
        # Directory case #
        if items[-1] == 'dir':
            items.pop(-1)
            return self.search_for_dir(key, items)
        else:
            return self.search_for_file(key, items)

    def search_for_file(self, key, items):
        # Search #
        matches = [set([p for p in self._paths if i in p]) for i in items]
        result = set.intersection(*matches)
        # No matches #
        if len(result) == 0:
            raise Exception("Could not find any path matching '%s'" % key)
        # Multiple matches, advantage file name #
        if len(result) > 1:
            best_score = max([p.score_file(items) for p in result])
            result = [p for p in result if p.score_file(items) >= best_score]
        # Multiple matches, take the one with less parts #
        if len(result) > 1:
            shortest = min([len(p) for p in result])
            result = [p for p in result if len(p) <= shortest]
        # Multiple matches, error #
        if len(result) > 1:
            raise Exception("Found several paths matching '%s'" % key)
        # Make the directory #
        result = result.pop()
        try:
            if not os.path.exists(result.complete_dir): os.makedirs(result.complete_dir)
        except OSError:
            pass
        # End #
        return result.complete_path

    def search_for_dir(self, key, items):
        # Search #
        matches = [set([p for p in self._paths if i in p]) for i in items]
        result = set.intersection(*matches)
        # No matches #
        if len(result) == 0:
            raise Exception("Could not find any path matching '%s'" % key)
        # Multiple matches, advantage dir name #
        if len(result) > 1:
            best_score = max([p.score_dir(items) for p in result])
            result = [p for p in result if p.score_dir(items) >= best_score]
        # Multiple matches, take the one with less parts #
        if len(result) > 1:
            shortest = min([len(p) for p in result])
            result = [p for p in result if len(p) <= shortest]
        # Multiple matches, maybe they all are the same directory #
        if len(result) > 1:
            if len(set([p.dir for p in result])) == 1: result = [result[0]]
        # Multiple matches, error #
        if len(result) > 1:
            raise Exception("Found several paths matching '%s'" % key)
        # Make the directory #
        result = result.pop()
        try:
            if not os.path.exists(result.complete_dir): os.makedirs(result.complete_dir)
        except OSError:
            pass
        # End #
        return result.complete_dir

    @property
    def tmp_dir(self):
        if not self._tmp_dir: self._tmp_dir = tempfile.mkdtemp() + '/'
        return self._tmp_dir

    @property
    def tmp(self):
        return self.tmp_dir + 'autopath.tmp'

################################################################################
class PathItems(object):
    delimiters = '_', '.', '/'
    pattern = re.compile('|'.join(map(re.escape, delimiters)))

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, path, base_dir):
        self.path = path
        self.base_dir = base_dir
        self.dir, self.name = os.path.split(path)
        self.name_items = self.pattern.split(self.name) if self.name else []
        self.dir_items = self.pattern.split(self.dir) if self.dir else []
        self.all_items = self.name_items + self.dir_items

    def __contains__(self, i):
        return i in self.all_items

    def __len__(self):
        return len(self.all_items)

    def score_file(self, items):
        return sum([1.0 if i in self.name_items else 0.5 for i in items if i in self])

    def score_dir(self, items):
        return sum([1.0 if i in self.dir_items else 0.5 for i in items if i in self])

    @property
    def complete_path(self):
        return '/' + os.path.relpath(self.base_dir + self.path, '/')

    @property
    def complete_dir(self):
        return '/' + os.path.relpath(self.base_dir + self.dir, '/') + '/'


################################################################################
class FilePath(object):

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, path):
        self.path = path

    @property
    def prefix_path(self):
        """The full path without the extension"""
        return os.path.splitext(self.path)[0]

    @property
    def prefix(self):
        """Just filename without the extension"""
        return os.path.basename(self.prefix_path)

    @property
    def directory(self):
        return os.path.dirname(self.path)

    @property
    def extension(self):
        return os.path.splitext(self.path)[1]
