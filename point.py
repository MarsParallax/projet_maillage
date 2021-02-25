"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np


class Point:
    """A simple point.

    Parameters
    ----------
    id : int
        a unique identifier
    X : numpy.array
        the coordinates
    tag : int
        the tag from GMSH
    """

    _counter = 1
    _name = "Point"

    def __init__(self, X, tag):
        """

        Parameters
        ----------
        X : numpy.array
            the coordinates
        tag : int
            the tag from GMSH
        """
        self.id = Point._counter
        Point._counter += 1
        self.X = X
        self.tag = tag

    # def __sub__(self, other):
    #     return Point(self.X - other.X)

    def distance(self, other):
        """Returns the distance between two Points.

        Parameters
        ----------
        other : Point
            the other point

        Returns
        -------
        float
            the distance
        """
        return np.linalg.norm(self.X - other.X)
