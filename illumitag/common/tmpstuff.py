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

    def __enter__(self):
        self.handle = open(self.path, 'w')
        return self

    def __exit__(self, *exc_info):
        self.handle.close()

    def __init__(self, path=None, **kwargs):
        if not path:
            handle = tempfile.NamedTemporaryFile(delete=False, **kwargs)
            self.path = handle.name
            handle.close()
        else:
            self.path = path

    def __repr__(self): return self.path
