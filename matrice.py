from triplets import Triplets

def assemblage(mesh):
    t = Triplets()
    for triangle in mesh.triangles:
        m_p = mat_elem(triangle)
        for i in range(2):
            I = l2g(triangle, i)
            for j in range(2):
                J = l2g(triangle, j)
                t.append(I, J, m_p(i, j))
    return t

def mat_elem(triangle)
    t = Triplets()
    area = triangle.area()
    for i in range(2):
        for j in range(2):
            if i == j:
                triplet.append(i, j, area / 6)
            else:
                triplet.append(i, j, area / 12)
    return t

def mass(mesh, dim, physical_tag, triplets): # TODO
    return triplets

def l2g(triangle, i):
    I = 0 #TODO
    return I