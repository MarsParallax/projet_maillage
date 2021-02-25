"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np

import point


class Triangle:
    """A simple triangle."""

    _counter = 1
    _name = "Triangle"

    def __init__(self, points, tag):
        self.id = Triangle._counter
        Triangle._counter += 1
        self.tag = tag
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

    def grad_phi_chap(self, i):
        grad = np.array([[-1, -1], [1, 0], [0, 1]])
        return grad[i]

    def matrice_rigidite_elem(self):
        D = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                BtB = np.dot(self.passage().T, self.passage())
                BtB_dot_grad = np.dot(self.grad_phi_chap(j).T, BtB)
                D[i][j] = self.area()*np.dot(BtB_dot_grad, self.grad_phi_chap(i))
        return D
