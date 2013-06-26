# Built-in modules #
import tempfile

################################################################################
class TmpFile(object):
    @classmethod
    def empty(cls, **kwargs):
        handle = tempfile.NamedTemporaryFile(delete=False, **kwargs)
        path = handle.name
        handle.close()
        return cls(path)

    @classmethod
    def from_string(cls, string, **kwargs):
        handle = tempfile.NamedTemporaryFile(delete=False, **kwargs)
        path = handle.name
        handle.write(string)
        handle.close()
        return cls(path)

    def __init__(self, path):
        self.path = path

    def __repr__(self): return self.path
