import point
import segment
import triangle
import mesh

class Assemblage :


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

	def loc2glob(node) :
		"""
			Giving a node 
		"""


