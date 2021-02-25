"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""

import numpy as np

import point


class Triangle:
    """Un simple triangle."""

    _counter = 1
    _name = "Triangle"

    def __init__(self, points, physical_tag=-1):
        self.id = Triangle._counter
        Triangle._counter += 1
        self.physical_tag = physical_tag
        assert len(points) == 3, "A triangle must have 3 points"
        self.points = points
        # self._area = 0.0
        self._set_area()

    def _set_area(self):
        x0 = self.points[0].X[0]
        y0 = self.points[0].X[1]
        x1 = self.points[1].X[0]
        y1 = self.points[1].X[1]
        x2 = self.points[2].X[0]
        y2 = self.points[2].X[1]
        self._area = abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0))

    def area(self):
        return self._area

    def jac(self):
        return 2 * self._area

    def passage(self):
        x1 = self.points[0].X[0]
        x2 = self.points[1].X[0]
        x3 = self.points[2].X[0]
        y1 = self.points[0].X[1]
        y2 = self.points[1].X[1]
        y3 = self.points[2].X[1]
        return (1 / self.jac()) * np.array([[y3 - y1, y1 - y2], [x1 - x3, x2 - x1]])
