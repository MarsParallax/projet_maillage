"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import sys

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse

from matrice import assemblage, dirichlet
from mesh import Mesh


def dirichlet_eval(X):
    """Evaluate Dirichlet condition on X.

    Parameters
    ----------
    X : numpy array
        the coordinates on which evaluation is done

    Returns
    -------
    double
        the result of the evaluation
    """
    return 0.0


# Load the mesh
# Init GMSH
gmsh.initialize(sys.argv)
# Ask GMSH to display information in the terminal
gmsh.option.setNumber("General.Terminal", 1)

filename = "mesh_project.msh"
model = gmsh.open(filename)
mesh = Mesh()
mesh.gmsh_to_mesh(model.mesh)

# Solve the problem
t, b = assemblage(mesh)
dirichlet(mesh, 1, 2, dirichlet_eval, t, b) # sur Rad
dirichlet(mesh, 1, 3, dirichlet_eval, t, b) # sur Fen
A = (scipy.sparse.coo_matrix(t.data)).tocsr()
U = scipy.sparse.linalg.spsolve(A, b)

# Plot the results
x = [point[0] for point in mesh.points]
y = [point[1] for point in mesh.points]
connectivity = []
for triangle in mesh.triangles:
    connectivity.append([point.tag for point in triangle.points])
plt.tricontourf(x, y, connectivity, U, 12)
plt.colorbar()
plt.show()

# Finalize GMSH
gmsh.finalize()
