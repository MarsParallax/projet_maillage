"""Copyright (c) Jérôme Bonacchi et Homer Durand 2021"""


from scipy.sparse import coo_matrix


class Triplets:
    """Un triplet de trois listes : valeurs, indices de lignes, indices de colonnes. 
    """

    def __init__(self):
        self.data = ([], ([], []))

    def __str__(self):
        return str(self.data)

    def append(self, I, J, val):
        """Ajoute le triplet [I, J, val] dans `self.data`.

        Parameters
        ----------
        I : [type]
            [description]
        J : [type]
            [description]
        val : [type]
            [description]
        """
        self.data[0].append(val)
        self.data[1][0].append(I)
        self.data[1][1].append(J)


if __name__ == "__main__":
    t = Triplets()
    t.append(0, 0, 1.1)
    t.append(0, 3, 2)
    t.append(1, 1, 1)
    t.append(2, 2, 2.3)
    t.append(3, 0, 0.5)
    t.append(3, 1, 2)
    t.append(3, 3, 2)
    A = coo_matrix(t.data)
    print(A.toarray())
    # [[1.1 0.  0.  2. ]
    # [0.  1.  0.  0. ]
    # [0.  0.  2.3 0. ]
    # [0.5 2.  0.  2. ]]
    A.tocsr()
