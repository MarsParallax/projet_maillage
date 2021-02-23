"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""

import gmsh

from point import Point
from segment import Segment
from triangle import Triangle


class Mesh:

    def gmsh_to_mesh(filename=None)
       if self.filename:
            print("Mesh already loaded")
        else:
            self.filename = filename
            self.N_points = 0  # TODO
            for (tag, X, _) in gmsh.model.mesh.getNodes():
                point = Point(X, tag)
                self.points.append(point)
            for (dim, physical_tag) in gmsh.model.mesh.getPhysicalGroups():
                for entity_tag in gmsh.model.mesh.getEntitiesForPhysicalGroup(dim, tag):
                    types, _, node_tags = gmsh.model.mesh.getElements(
                        dim, entity_tag)
                    print(types, node_tags)
                    # types = gmsh.model.mesh.getElementTypes(dim,entity)
                    # if types[0] == 1:
                    #     X =
                    #     segment = Segment(X, tag)
                    #     self.segments.append(segment)
                    # elif types[0] == 2:
                    #     X =
                    #     triangle = Triangle(X, tag)
                    #     self.triangles.append(triangle)
                    # else:
                    #     print("Unknown type.")

    def get_points(self, dim, physical_tag):
        """Return: nodeTags, coord"""
        gmsh.model.mesh.getNodesForPhysicalGroup(dim, physical_tag)

    def get_elements(self, dim, physical_tag):  # TODO
        tags = gmsh.model.mesh.getEntitiesForPhysicalGroup(dim, physical_tag)
        # elements = list()
        # for tag in tags:
        #     gmsh.model.mesh.getElements(dim, tag)
        return tags
