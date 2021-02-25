"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np


class Point:
    """A simple point."""

    _counter = 1
    _name = "Point"

    def __init__(self, X, tag):
        self.id = Point._counter
        Point._counter += 1
        self.X = X
        self.tag = tag

    # def __sub__(self, other):
    #     return Point(self.X - other.X)

    def distance(self, other):
        return np.linalg.norm(self.X - other.X)
