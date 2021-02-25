import point
import segment
import triangle
import mesh
import numpy as np


def assembler() :
	A = 0
	b = 0
	For p = 1:N_t
	  For i = 1:3
	    I = locToGlob(p,i)
	    For j = 1:3
	      J = locToGlob(p,j)
	      A(I,J) += ∫_{K_p}(∇ϕ_j^p·∇ϕ_i^p)
	    EndFor
	    b(I) += ∫_{K_p}(∇u_ɣ^h·∇ϕ_i^p)
	  EndFor
	EndFor
	return 0



def rigidite_elem(element) :
	return 0

def loc2glob(p, i) :
	"""
		 
	"""

def grad_phi_chap(p, i) :
	"""

	"""
	grad = np.array([[-1, -1], [1, 0], [0, 1]])
	return grad_phi[i]
