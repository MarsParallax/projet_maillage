"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import sys

import gmsh
import numpy as np
from matplotlib import pyplot as plt

from point import Point
from segment import Segment
from triangle import Triangle


class Mesh:
    """Wraps the calls to GMSH"""

    def __init__(self):
        self.points = []
        self.segments = []
        self.triangles = []
        self.filename = ""

    def gmsh_to_mesh(self, filename):
        """Initialize self.points, self.segments and self.triangles with mesh
        stored in filename

        Parameters
        ----------
        filename : str
            the file name
        """
        self.filename = filename
        gmsh.open(filename)

        # Fill self.points
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        X = nodeCoords[0::3]
        Y = nodeCoords[1::3]
        for (tag, x, y) in zip(nodeTags, X, Y):
            point = Point(np.array([x, y]), tag)
            self.points.append(point)

        # Set number of points of mesh
        self.nbPoints = len(X)

        # Fill self.segments, self.triangles
        for dim, physical_tag in gmsh.model.getPhysicalGroups():
            for entity_tag in gmsh.model.getEntitiesForPhysicalGroup(dim, physical_tag):
                elmType, elmTags, nodeTags = gmsh.model.mesh.getElements(
                    dim, entity_tag)

                # Segment
                if dim == 1:
                    for j in range(len(elmTags[0])):
                        ids = [nodeTags[0][2 * j], nodeTags[0][2 * j + 1]]
                        pts = [point for point in self.points if point.id in ids]
                        segment = Segment(pts, elmTags[0][j])
                        self.segments.append(segment)

                # Triangle
                elif dim == 2:
                    for j in range(len(elmTags[0])):
                        ids = [nodeTags[0][3 * j], nodeTags[0][3 * j + 1],
                               nodeTags[0][3 * j + 2]]
                        points = [
                            point for point in self.points if point.tag in ids]
                        try:
                            triangle = Triangle(points, elmTags[0][j])
                            self.triangles.append(triangle)
                        except Exception as e:
                            print("[WARNING] `Mesh.gmsh_to_mesh`:", e)
                            print("Points are:", [point.X for point in points])
                # Others
                else:
                    print(
                        f"[WARNING] `Mesh.gmsh_to_mesh`: Unknown element type '{dim}'.")

    def get_points(self, dim=-1, physical_tag=-1):
        """Return points in the entity of dimensions dim et physical tag
        physical_tag. By default, returns all the points.

        Parameters
        ----------
        dim : int
            the dimension, by default -1
        physical_tag : int
            the physical_tag, by default -1

        Returns
        -------
        numpy.array
            the points
        """
        if dim == - 1 and physical_tag == - 1:
            return np.array([self.points])
        # Return points contained in segments
        if dim == 1:
            if physical_tag == - 1:
                ids = np.array(
                    [point.id for segment in self.segments for point in segment.points])
            else:
                ids = np.array(
                    [point for segment in self.segments for point in segment.points if segment.tag == physical_tag])
            return np.array([point for point in self.points if point.id in np.unique(ids)])
        # Returns points contained in triangles
        elif dim == 2:
            if physical_tag == - 1:
                ids = np.array(
                    [point.id for triangle in self.triangles for point in triangle.points])
            else:
                ids = np.array(
                    [point.id for triangle in self.triangles for point in triangle.points if triangle.tag == physical_tag])
            return np.array([point for point in self.points if point.id in np.unique(ids)])
        # Others
        else:
            print(
                f"[WARNING] `Mesh.get_points`: Unknown element type '{dim}'.")
        return None

    def get_elements(self, dim, physical_tag):
        """Return elements of dimension dim and physical tag physical_tag
        (dim = 1 for segement and 2 for triangles).

        Parameters
        ----------
        dim : int
            the dimension
        physical_tag : int
            the physical tag

        Returns
        -------
        numpy.array
            the elements
        """
        # Return segments
        if dim == 1:
            if physical_tag == -1:
                return np.array([segment for segment in self.segments])
            else:
                return np.array([segment for segment in self.segments if segment.tag == physical_tag])
        # Return triangles
        elif dim == 2:
            if physical_tag == -1:
                return np.array([triangle for triangle in self.triangles])
            else:
                return np.array([triangle for triangle in self.triangles if triangle.tag == physical_tag])
        # Others
        else:
            print(
                f"[WARNING] `Mesh.get_elements`: Unknown element type '{dim}'.")
        return None

    def __str__(self):
        res = [f"Mesh loaded from '{self.filename}' composed of:",
               f"  {len(self.points)} points",
               f"  {len(self.segments)} segments",
               f"  {len(self.triangles)} triangles"]
        return "\n".join(res)


if __name__ == '__main__':
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    mesh = Mesh()
    mesh.gmsh_to_mesh("domaine.msh")
    print(mesh)
    elm = mesh.get_elements(1, 2)
    pt = mesh.get_points(2, 2)
    gmsh.finalize()
