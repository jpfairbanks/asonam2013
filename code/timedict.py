from time import time as time

class timedict(object):
    """ The global state necessary to do sophisticated timing
    acts like a key value store of times using tic and toc to start and stop the timers.

    timedict.as_dict() returns the dict from keys to the time they took. the
    keys can be any hashable type. Such as a loop iteration counter or function
    name as a string
    """

    def __init__(self, ):
        """
        """
        self.starts = dict()
        self.ends = dict()

    def tic(self,  key):
        """start the timer with key

        Arguments:
        - `self`: state
        - `key`: the key to access this later
        """
        self.starts[key] = time()

    def toc(self, key):
        """start the timer with key

        Arguments:
        - `self`: state
        - `key`: the key to access this later
        """
        self.ends[key] = time()-self.starts[key]

    def __str__(self):
        return self.ends.__str__()

    def as_dict(self):
        return self.ends
