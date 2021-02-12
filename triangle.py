import point

class Triangle:

    _counter = 0
    _name = "Triangle"

    def __init__(self, points, physical_tag=-1)
        self.id = Trangle._counter
        Triangle_counter += 1
        self.physical_tag = physical_tag
        assert len(points) == 3, "A triangle must have 3 points"
        self.points = points
        self._area = 0.0
        self._set_area()

    def _set_area(self):
        x0 = self.points[0][0]
        y0 = self.points[0][1]
        x1 = self.points[1][0]
        y1 = self.points[1][1]
        x2 = self.points[2][0]
        y2 = self.points[2][1]
        self._area = abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0))

    def area(self):
        return self._area

    def jac(self):
        return 2 * self._area
