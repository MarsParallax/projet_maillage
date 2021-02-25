"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np

import point


class Segment:
    """A simple segment."""

    _counter = 1
    _name = "Segment"

    def __init__(self, points, tag):
        """Initialize a sement with a set of points. 
        Parameters
        ----------
        points : list of Points
            points that represent the segement
        tag : int
            tag of the point
        """
        self.id = Segment._counter
        Segment._counter += 1
        self.tag = tag
        assert len(points) == 2, "A segment must have 2 points"
        self.points = points
        # self._length = 0.0
        self._set_length()

    def _set_length(self):
        """
            Initialize the length of the segment
        """
        self._length = self.points[0].distance(self.points[1])

    def area(self):
        """
            Return the length of the segment 
        """
        return self._length

    def jac(self):
        """Return jacobian 
        """
        return self._length
