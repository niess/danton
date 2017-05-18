# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
# Author: Valentin NIESS (niess@in2p3.fr)
#
# This software is a C99 executable dedicated to the sampling of decaying
# taus from ultra high energy neutrinos.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

import collections

# The Earth radius in the Preliminary Earth Model (PEM).
EARTH_RADIUS = 6371.E+03

# Data structures for decays.
State = collections.namedtuple("State", ("pid", "energy", "direction",
    "position"))
Product = collections.namedtuple("Product", ("pid", "momentum"))
Decay = collections.namedtuple("Decay", ("generation", "tau_i", "tau_f",
    "product"))
Event = collections.namedtuple("Event", ("id", "primary", "decay", "weight"))

class iter_event:
    """Iterator over the tau decays in a text dump.
    """
    def __init__(self, filename):
        self.fid = open(filename, "r")
        for _ in xrange(3): self.fid.readline() # skip the header
        self.field = self.fid.readline().split()

    def __iter__(self):
        return self

    def next(self):
        if not self.field:
                if not self.fid:
                        self.fid.close()
                        self.fid = None
                raise StopIteration()
        self.field, event = self._get_next_event()
        return event

    def _get_next_event(self):
        """Get the next event in the record.
        """
        # Get the primary and event info.
        eventid = int(self.field[0])
        weight = float(self.field[9])
        primary = State(int(self.field[1]),
            float(self.field[2]), map(float, self.field[3:6]),
            map(float, self.field[6:9]))

        # Get the tau decays info.
        decay = []
        self.field = self.fid.readline().split()
        while len(self.field) == 9:
                genid = int(self.field[0])
                pid = int(self.field[1])
                tau_i = State(pid, float(self.field[2]),
                    map(float, self.field[3:6]), map(float, self.field[6:9]))
                self.field = self.fid.readline().split()
                tau_f = State(pid, float(self.field[0]),
                    map(float, self.field[1:4]), map(float, self.field[4:7]))

                # Get the decay product(s)
                product = []
                self.field = self.fid.readline().split()
                while len(self.field) == 4:
                        product.append(Product(int(self.field[0]),
                                        map(float, self.field[1:4])))
                        self.field = self.fid.readline().split()
                decay.append(Decay(genid, tau_i, tau_f, product))
        return self.field, Event(eventid, primary, decay, weight)
