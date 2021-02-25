"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


from itertools import product

import numpy as np
from scipy.sparse import coo_matrix

import mesh
from triplets import Triplets


def assemblage(mesh):
    """Assembles A and b

    Parameters
    ----------
    mesh : Mesh
        the mesh

    Returns
    -------
    (Triplets, numpy.array)
        matrix A and vector b
    """
    A = Triplets()
    b = np.zeros((mesh.nbPoints))
    for p in mesh.get_elements(2, -1):
        # for triangle in mesh.triangles:
        Dp = p.matrice_rigidite_elem()
        bp = p.rhs()
        for i in range(3):
            I = local_to_global(p, i)
            for j in range(3):
                J = local_to_global(p, j)
                A.append(I, J, Dp[i][j])
            b[I] += bp[i]
    return (A, b)

# def mass_elem(triangle):
#     """Computes the elementary mass matrix on a triangle

#     Parameters
#     ----------
#     triangle : Triangle
#         the triangle on which the matrix is computed

#     Returns
#     -------
#     Triplets (length = 9)
#         the elementary mass matrix
#     """
#     M = Triplets()
#     area = triangle.area()
#     for (i, j) in product(range(2), repeat=2):
#         if i == j:
#             M.append(i, j, area / 6)
#         else:
#             M.append(i, j, area / 12)
#     return M


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
    I = triangle.points[i].id - 1
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
    corresp = {int(point.id - 1): point for point in points}
    I = [int(point.id - 1) for point in points]
    row_indices = triplets.data[1][0]
    for i in row_indices:
        if i in I:
            triplets.data[0][i] = 0
    for i in I:
        triplets.append(i, i, 1)
        b[i] = g(corresp[i].X)


if __name__ == '__main__':
    mesh = mesh.Mesh()
    mesh.gmsh_to_mesh("domaine.msh")
    A, _ = assemblage(mesh)
    A_coo = coo_matrix(A.data).tocsr()
    U = np.zeros((A_coo.get_shape()[0]))+1
    print(A_coo.dot(U))
