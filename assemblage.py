import point
import segment
import triangle
import mesh

class Assemblage :


	def assembler() :
		Triplets t;
		For p = 1, ... Nt  // Parcours des Triangles
		  Mp = MatElem(p); // Matrice Elementaire du triangle p
		  For i = 1,2,3
		    I = Loc2Glob(p, i);
		    For j = 1,2,3
		      J = Loc2Glob(p,j);
		      t.append(I, J, Mp(i,j)); // contribution élémentaire
		    End
		End
		return 0

	def rigidite_elem(element) :
		return 0

	def loc2glob() :


