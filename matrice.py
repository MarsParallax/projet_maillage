"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


from itertools import product
from triplets import Triplets
import mesh
import numpy as np
from scipy.sparse import coo_matrix

def assemblage(mesh) :
	""" Assemble A

    Parameters
    ----------
    mesh : Mesh
        the mesh

    Returns
    -------
    Triplets
        matrix A
    numpy.array
    	matrix B
    """
	A = Triplets()
	b = np.zeros((mesh.nbPoints))
	for p in mesh.get_elements(2, -1) :
		Mp = p.matrice_rigidite_elem()
		for i in range(3) :
			I = local_to_global(p,i)
			for j in range(3) :
				J = local_to_global(p,j)
				A.append(I, J, Mp[i][j])
			#b[I] += 0
	return A, b



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


if __name__ == '__main__':
    mesh = mesh.Mesh()
    mesh.gmsh_to_mesh("domaine.msh")
    A, _ = assemblage(mesh) 
    A_coo = coo_matrix(A.data).tocsr()
    U = np.zeros((A_coo.get_shape()[0]))+1
    print(A_coo.dot(U))