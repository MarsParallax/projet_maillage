import gmsh
import sys

import matplotlib.pyplot as plt
import numpy as np

from mesh import Mesh
from matrice import *

# Init GMSH
gmsh.initialize(sys.argv)
# Ask GMSH to display information in the terminal
gmsh.option.setNumber("General.Terminal", 1)

filename = "mesh_project.msh"
model = gmsh.open(filename)
mesh = Mesh()
mesh.gmsh_to_mesh(model.mesh)
t = Triplets()
mass(mesh, , , t)
b = np.zeros((mesh.Npts,))
A = (sparse.coo_matrix(t.data)).tocsr()
U = sparse.linalg.spsolve(A, b)

x = [point[0] for point in mesh.points]
y = [point[1] for point in mesh.points]
connectivity=[]
for triangle in mesh.triangles:
    connectivity.append([point.id for point in triangle.points])
plt.tricontourf(x, y, connectivity, U, 12)
plt.colorbar()
plt.show()

# Finalize GMSH
gmsh.finalize()