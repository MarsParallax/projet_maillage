"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


from itertools import product

import numpy as np

from triplets import Triplets


def assemblage(mesh):
    """Assembles the global matrix.

    Parameters
    ----------
    mesh : Mesh
        the mesh

    Returns
    -------
    Triplets
        the matrix as a Triplets object
    """
    A = Triplets()
    b = np.zeros((mesh.Npts,))
    for triangle in mesh.triangles:
        # Mp = mass_elem(triangle)
        Dp = stiffness_elem(triangle)
        for i in range(2):
            I = local_to_global(triangle, i)
            for j in range(2):
                J = local_to_global(triangle, j)
                A.append(I, J, Dp[i, j])
            b[I] += rhs(triangle, i)
    return (A, b)


def rhs(triangle, i):
    """Computes the right-hand side

    Parameters
    ----------
    triangle : Triangle
        the triangle
    i : int
        the local index

    Returns
    -------
    float
        the value
    """
    return triangle.area # TODO


def mass_elem(triangle):
    """Computes the elementary mass matrix on a triangle

    Parameters
    ----------
    triangle : Triangle
        the triangle on which the matrix is computed

    Returns
    -------
    Triplets (length = 9)
        the elementary mass matrix
    """
    M = Triplets()
    area = triangle.area()
    for (i, j) in product(range(2), repeat=2):
        if i == j:
            M.append(i, j, area / 6)
        else:
            M.append(i, j, area / 12)
    return M


def stiffness_elem(triangle):
    """Computes the elementary stiffness matrix on a triangle

    Parameters
    ----------
    triangle : Triangle
        the triangle on which the matrix is computed

    Returns
    -------
    Triplets (length = 9)
        the elementary stiffness matrix
    """
    D = Triplets()
    # area = triangle.area()
    for (i, j) in product(range(2), repeat=2):
        if i == j:
            D.append(i, j, 0)  # TODO
        else:
            D.append(i, j, 0)  # TODO
    return D


def local_to_global(triangle, i):
    """Converts the local index into the global index.

    Parameters
    ----------
    triangle : Triangle
        the triangle
    i : int
        the local index

    Returns
    -------
    int
        the global index
    """
    I = 3 * triangle.id + i
    return I


def dirichlet(mesh, dim, physical_tag, g, triplets, b):
    """Applies the Dirichlet condition u=g on the domain corresponding to dim
    and physical tag in the system made of triplet and b. triplets and b are
    altered.

    Parameters
    ----------
    mesh : Mesh
        the mesh
    dim : int
        the dimension of the domaine on which the condition is applied
    physical_tag : int
        the physical tag of the domaine on which the condition is applied
    g : function(numpy.array)
        the function of the Dirichlet condition
    triplets : Triplets
        the matrix of the system
    b : numpy.array
        the right-hand side of the system
    """
    points = mesh.get_points(dim, physical_tag)
    I = [point.tag for point in points]
    values =  triplets.data[0]
    row_indices = triplets.data[1][0]
    for i in row_indices: 
        if i in I:
            values[i] = 0
    for i in I:
        triplets.append(i, i, 1)
        b[i] = g(points[i].X)
