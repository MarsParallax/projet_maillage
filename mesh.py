"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import sys

import gmsh
import numpy as np 
from matplotlib import pyplot as plt

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

from point import Point
from segment import Segment
from triangle import Triangle


class Mesh:
    """Wraps the gmsh calls"""

    def __init__(self):
        self.points = []
        self.segments = []
        self.triangles = []

    def gmsh_to_mesh(self, filename):
        """Initialize self.points, self.segments and self.triangles with mesh
        stored in filename

        Parameters
        ----------
        filename : str
            the file name
        """
        gmsh.open(filename)

        # Fill self.points
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        X = nodeCoords[0::3]
        Y = nodeCoords[1::3]
        for (tag, x, y) in zip(nodeTags, X, Y):
            point = Point(np.array([x, y]), tag)
            self.points.append(point)

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
                        ids = [nodeTags[0][2 * j], nodeTags[0][2* j + 1],
                               nodeTags[0][2 * j + 2]]
                        X = [point for point in self.points if point.id in ids]
                        try:  # TODO : pb triangle plat
                            triangle = Triangle(X, elmTags[0][j])
                            self.triangles.append(triangle)
                        except Exception as e:
                            print("[WARNING] `Mesh.gmsh_to_mesh`:", e)
                # Others
                else:
                   print(f"[WARNING] `Mesh.gmsh_to_mesh`: Unknown element type '{dim}'.")

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
                ids = np.array([point.id for segment in self.segments for point in segment.points])
            else:
                ids = np.array([point for segment in self.segments for point in segment.points if segment.physical_tag == physical_tag])
            return np.array([point for point in self.points if point.id in np.unique(ids)])
        # Returns points contained in triangles
        elif dim == 2:
            if physical_tag == - 1:
                ids = np.array([point.id for triangle in self.triangles for point in triangle.points])
            else:
                ids = np.array([point.id for triangle in self.triangles for point in triangle.points if triangle.physical_tag == physical_tag])
            return np.array([point for point in self.points if point.id in np.unique(ids)])
        # Others
        else:
            print(f"[WARNING] `Mesh.get_points`: Unknown element type '{dim}'.")
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
                return np.array([segment for segment in self.segments if segment.physical_tag == physical_tag])
        # Return triangles
        elif dim == 2:
            if physical_tag == -1:
                return np.array([triangle for triangle in self.triangles])
            else:
                return np.array([triangle for triangle in self.triangles if triangle.physical_tag == physical_tag])
        # Others
        else:
            print(f"[WARNING] `Mesh.get_elements`: Unknown element type '{dim}'.")
        return None


if __name__ == '__main__':
    mesh = Mesh()
    mesh.gmsh_to_mesh("domaine.msh")
    elm = mesh.get_elements(1, 2)
    pt = mesh.get_points(2, 2)
