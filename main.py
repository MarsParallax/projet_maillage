"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import sys

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

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

filename = "domaine.msh"
mesh = Mesh()
try:
    mesh.gmsh_to_mesh("domaine.msh")
except AttributeError as e:
    print("[WARNING] `main.py`: check if the file name is correct.")
    print(e)

print("nb triangle : ", len(mesh.get_elements(2, -1)))
# Solve the problem
t, b = assemblage(mesh)
dirichlet(mesh, 1, 2, dirichlet_eval, t, b) # sur Rad
dirichlet(mesh, 1, 3, dirichlet_eval, t, b) # sur Fen
A = (scipy.sparse.coo_matrix(t.data)).tocsr()
U = spsolve(A, b)

# Plot the results
x = [point.X[0] for point in mesh.points]
y = [point.X[1] for point in mesh.points]
connectivity = []


for triangle in mesh.triangles:
    connectivity.append([int(point.tag-1) for point in triangle.points])

fig, ax = plt.subplots()  # figsize=(12, 6)
ax.tricontour(x, y, connectivity, U, levels=12, linewidths=0.5, colors='k')
contour = ax.tricontourf(x, y, connectivity, U, levels=12, cmap="RdBu_r")
fig.colorbar(contour, ax=ax)
ax.axis("scaled")
ax.set(xlim=(0, 10), ylim=(0, 10), title="Simulation")
plt.show()

# Finalize GMSH
gmsh.finalize()