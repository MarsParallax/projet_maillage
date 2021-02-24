"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np


class Point:
    """Un simple point.

    Returns
    -------
    [type]
        [description]
    """

    _counter = 1
    _name = "Point"

    def __init__(self, X, global_index=None):
        self.id = Point._counter
        Point._counter += 1
        self.X = X

    def __sub__(self, other):
        return Point(self.X - other.X)

    def distance(self, other):
        return np.linalg.norm(self.X - other.X)
