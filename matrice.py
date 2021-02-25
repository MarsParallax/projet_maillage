"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


from itertools import product

from triplets import Triplets


def assemblage(mesh):
    """Assemble la matrice du système.

    Parameters
    ----------
    mesh : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    A = Triplets()
    for triangle in mesh.triangles:
        # Mp = mass_elem(triangle)
        Dp = rigid_elem(triangle)
        for i in range(2):
            I = local_to_global(triangle, i)
            for j in range(2):
                J = local_to_global(triangle, j)
                A.append(I, J, Dp[i, j])
    return A


def mass_elem(triangle):
    """Calcule la matrice de masse élémentaire d'un triangle

    Parameters
    ----------
    triangle : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    M = Triplets()
    area = triangle.area()
    for (i, j) in product(range(2), repeat=2):
        if i == j:
            M.append(i, j, area / 6)
        else:
            M.append(i, j, area / 12)
    return M

def rigid_elem(triangle):
    """Calcule la matrice de rigidité élémentaire d'un triangle

    Parameters
    ----------
    triangle : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    M = Triplets()
    area = triangle.area()
    for (i, j) in product(range(2), repeat=2):
        if i == j:
            M.append(i, j, 0) # TODO
        else:
            M.append(i, j, 0) # TODO
    return M


def mass(mesh, dim, physical_tag, triplets):  # TODO
    return triplets


def local_to_global(triangle, i):
    """Convertit l'indice local à un triangle en l'indice global.

    Parameters
    ----------
    triangle : [type]
        [description]
    i : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    I = 3 * triangle.physical_tag + i # TODO
    return I
