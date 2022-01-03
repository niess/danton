'''Python utilities for the Danton Library
'''
import collections
try:
    from cStringIO import StringIO as BytesIO
except ImportError:
    # Python 3
    from io import BytesIO
else:
    # Python 2 compatibility
    range = xrange


# The Earth radius in the Preliminary Earth Model (PEM).
EARTH_RADIUS = 6371.E+03


class iter_event:
    '''Iterator over the tau decays in a text dump.
    '''

    # Data structures for decays.
    State = collections.namedtuple('State', ('pid', 'energy', 'direction',
        'position'))
    Product = collections.namedtuple('Product', ('pid', 'momentum'))
    Decay = collections.namedtuple('Decay', ('generation', 'tau_i', 'tau_f',
        'product'))
    Event = collections.namedtuple('Event', ('id', 'primary', 'decay',
        'weight'))

    def __init__(self, filename=None, text=None):
        if text is not None:
            self.fid = BytesIO(text)
        else:
            self.fid = open(filename, 'r')
        for _ in range(3): self.fid.readline() # skip the header
        self.field = self.fid.readline().split()

    def __delete__(self):
        if self.fid is not None:
            self.fid.close()
            self.fid = None

    def __iter__(self):
        return self

    def next(self):
        if not self.field:
            if not self.fid:
                    self.fid.close()
                    self.fid = None
            raise StopIteration()
        try:
            self.field, event = self._get_next_event()
        except:
            raise StopIteration()
        return event

    def _get_next_event(self):
        '''Get the next event in the record.
        '''
        # Get the primary and event info.
        eventid = int(self.field[0])
        weight = float(self.field[9])
        primary = self.State(int(self.field[1]),
            float(self.field[2]), map(float, self.field[3:6]),
            map(float, self.field[6:9]))

        # Get the tau decays info.
        decay = []
        self.field = self.fid.readline().split()
        while len(self.field) == 9:
                genid = int(self.field[0])
                pid = int(self.field[1])
                tau_i = self.State(pid, float(self.field[2]),
                    map(float, self.field[3:6]), map(float, self.field[6:9]))
                self.field = self.fid.readline().split()
                tau_f = self.State(pid, float(self.field[0]),
                    map(float, self.field[1:4]), map(float, self.field[4:7]))

                # Get the decay product(s)
                product = []
                self.field = self.fid.readline().split()
                while len(self.field) == 4:
                        product.append(self.Product(int(self.field[0]),
                                        map(float, self.field[1:4])))
                        self.field = self.fid.readline().split()
                decay.append(self.Decay(genid, tau_i, tau_f, product))
        return self.field, self.Event(eventid, primary, decay, weight)


class iter_flux:
    '''Iterator over the flux events in a text dump.
    '''

    # Data structures for events.
    State = collections.namedtuple('State', ('pid', 'generation', 'energy',
        'direction', 'position'))
    Event = collections.namedtuple('Event', ('id', 'primary', 'final',
        'weight'))

    def __init__(self, filename=None, text=None):
        if text is not None:
            self.fid = BytesIO(text)
        else:
            self.fid = open(filename, 'r')
        for _ in range(3): self.fid.readline() # skip the header
        self.field = self.fid.readline().split()

    def __delete__(self):
        if self.fid is not None:
            self.fid.close()
            self.fid = None

    def __iter__(self):
        return self

    def __next__(self):
        if not self.field:
            if not self.fid:
                    self.fid.close()
                    self.fid = None
            raise StopIteration()
        self.field, event = self._get_next_event()
        return event

    def next(self):
        return self.__next__()

    def _get_next_event(self):
        '''Get the next event in the record.
        '''
        # Get the primary and event info.
        eventid = int(self.field[0])
        weight = float(self.field[9])
        primary = self.State(int(self.field[1]), 1,
            float(self.field[2]), map(float, self.field[3:6]),
            map(float, self.field[6:9]))

        # Get the final state(s) info.
        final = []
        self.field = self.fid.readline().split()
        if len(self.field) == 9:
            while len(self.field) == 9:
                genid = int(self.field[0])
                pid = int(self.field[1])
                if abs(pid) == 15:
                    self.field = self.fid.readline().split()
                    final.append(self.State(pid, genid, float(self.field[0]),
                        map(float, self.field[1:4]),
                        map(float, self.field[4:7])))
                else:
                    final.append(self.State(pid, genid, float(self.field[2]),
                        map(float, self.field[3:6]),
                        map(float, self.field[6:9])))

                self.field = self.fid.readline().split()

        return self.field, self.Event(eventid, primary, final, weight)
