"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


import sys

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from matrice import assemblage, dirichlet
from mesh import Mesh


def dirichlet_fen(X): return -10
def dirichlet_rad(X): return 25


# Load the mesh
# Init GMSH
gmsh.initialize(sys.argv)
# Ask GMSH to display information in the terminal
gmsh.option.setNumber("General.Terminal", 1)

argv = sys.argv
argc = len(argv)
if argc == 1:
    filename = "domaine_h1.msh"
elif argc == 2:
    filename = argv[1]
elif argc > 2:
    sys.exit("Usage: [filename]")
mesh = Mesh()
try:
    mesh.gmsh_to_mesh(filename)
    # Solve the problem
    t, b = assemblage(mesh)
    dirichlet(mesh, 1, 2, dirichlet_rad, t, b)  # sur Rad
    dirichlet(mesh, 1, 3, dirichlet_fen, t, b)  # sur Fen
    np.set_printoptions(threshold=sys.maxsize)
    A = (coo_matrix(t.data)).tocsr()
    U = spsolve(A, b)
except Exception as e:
    print("[WARNING] `main.py`: check first if the file name is correct.")
    print(e)

# Plot the results
x = [point.X[0] for point in mesh.points]
y = [point.X[1] for point in mesh.points]
connectivity = []
for triangle in mesh.triangles:
    connectivity.append([int(point.id-1) for point in triangle.points])

fig, ax = plt.subplots(figsize=(7, 7))
ax.tricontour(x, y, connectivity, U, levels=12, linewidths=0.5, colors='k')
contour = ax.tricontourf(x, y, connectivity, U, levels=12, cmap="RdBu_r")
fig.colorbar(contour, ax=ax)
ax.axis("scaled")
ax.set(xlim=(0, 10), ylim=(0, 10),
       title="Simulation de la température (en degré Celsius)")  # \n h=0,1
plt.show()

# Finalize GMSH
gmsh.finalize()
