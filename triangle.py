"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import numpy as np

import point


class Triangle:
    """A simple triangle.

    Parameters
    ----------
    id : int
        a unique identifier
    points : numpy.array
        the points of the triangle
    tag : int
        the tag from GMSH
    """

    _counter = 1
    _name = "Triangle"

    def __init__(self, points, tag):
        """Initializes Triangle with a set of 3 points and a tag

        Parameters
        ----------
        points : numpy.array
            the 3 points of the triangle
        tag : int
            the tag of the triangle from GMSH
        """
        self.id = Triangle._counter
        Triangle._counter += 1
        self.tag = tag
        assert len(points) == 3, "A triangle must have 3 points"
        self.points = points
        # self._area = 0.0
        self._set_area()

    def _set_area(self):
        """Initialize self.area, the area of the triangle
        """
        x0 = self.points[0].X[0]
        y0 = self.points[0].X[1]
        x1 = self.points[1].X[0]
        y1 = self.points[1].X[1]
        x2 = self.points[2].X[0]
        y2 = self.points[2].X[1]
        self._area = abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0))

    def area(self):
        """Return the area of the triangle 
        """
        return self._area

    def jac(self):
        """Return the jacobian of the triangle
        """
        return 2 * self._area

    def passage(self):
        """Return passage matrix from local gradient of phi to global gradient of phi
        """
        x1 = self.points[0].X[0]
        x2 = self.points[1].X[0]
        x3 = self.points[2].X[0]
        y1 = self.points[0].X[1]
        y2 = self.points[1].X[1]
        y3 = self.points[2].X[1]
        return (1 / self.jac()) * np.array([[y3 - y1, y1 - y2], [x1 - x3, x2 - x1]])

    def grad_phi_chap(self, i):
        """Return local gradient of phi.

        Parameters
        ----------
        i : int
            local id of the point
        """
        grad = np.array([[-1, -1], [1, 0], [0, 1]])
        return grad[i]

    def matrice_rigidite_elem(self):
        """return elementary rigidity matrix
            
        """
        D = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                BtB = np.dot(self.passage().T, self.passage())
                BtB_dot_grad = np.dot(self.grad_phi_chap(j).T, BtB)
                D[i][j] = self.area()*np.dot(BtB_dot_grad, self.grad_phi_chap(i))
        return D

    def rhs(self):
        b = np.zeros((3))
        for i in range(3):
            BtB = np.dot(self.passage().T, self.passage())
            BtB_dot_grad = np.dot(np.array([25,25]).T, BtB)
            b[i] = self.area()*np.dot(BtB_dot_grad, self.grad_phi_chap(i))
        return b
