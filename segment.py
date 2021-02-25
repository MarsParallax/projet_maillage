"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np

import point


class Segment:
    """A simple segment.
    
    Parameters
    ----------
    id : int
        a unique identifier
    tag : int
        the tag from GMSH
    points : numpy.array(Points)
        the two points of the segment
    length : float
        length of the segment
    """

    _counter = 1
    _name = "Segment"

    def __init__(self, points, tag):
        """Initializes a segment with a 2 points.

        Parameters
        ----------
        points : numpy.array(Points)
            the two points of the segment
        tag : int
            the tag from GMSH
        """
        self.id = Segment._counter
        Segment._counter += 1
        self.tag = tag
        assert len(points) == 2, "A segment must have 2 points"
        self.points = points
        # self._length = 0.0
        self._set_length()

    def _set_length(self):
        """Initializes the length of the segment.
        """
        self._length = self.points[0].distance(self.points[1])

    def area(self):
        """Returns the length of the segment.
        """
        return self._length

    def jac(self):
        """Returns jacobian.
        """
        return self._length
