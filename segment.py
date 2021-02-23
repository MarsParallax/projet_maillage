"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np

import point


class Segment:
    """Un simple segment."""

    _counter = 0
    _name = "Segment"

    def __init__(self, points, physical_tag=-1):
        self.id = Segment._counter
        Segment._counter += 1
        self.physical_tag = physical_tag
        assert len(points) == 2, "A 'segment must have 2 points"
        self.points = points
        # self._length = 0.0
        self._set_length()

    def _set_length(self):
        self._length = self.points[0].distance(self.points[1])

    def area(self):
        return self._length

    def jac(self):
        return self._length
