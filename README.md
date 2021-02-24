# Projet de maillage et éléments finis

par Jérôme Bonacchi et Homer Durand à Polytech Sorbonne en spécialité mathématiques appliquées et informatique

## TODO

- 1a
  - [ ] gmshToMesh : construire la listes des éléments
  - [ ] getElements, getPoints : faire un find, pour retourner les objets
- 1b
  - [x] modélisation en gmsh
- 1c
  - [ ] récrire le problème/la reformulation au propre
  - [ ] faire le relèvement
- 2
  - [ ] tester la modélisation gmsh
  - [ ] matrice de rigidité : calculer les gradPhi et la jacobienne https://bthierry.pages.math.cnrs.fr/course-fem/lecture/elements-finis-triangulaires/contributions-elementaires/
  - [ ] assemblage des matrices élémentaires
- 3
  - [ ] faire `main.py`
  - [ ] affichage graphique avec gradient de couleur
- 4
  - [ ] commenter le code
  - [ ] faire le readme : définir le bon Vh, vérifier la formulation

## Exécution

Le code a été testé avec :

- python 3.9.1
- scipy 1.6.0
- numpy 1.20
- matplotlib 3.3.4
- gmsh 4.7.1.

Pour lancer le programme, il suffit de tapper :

```sh
python3 main.py
```

## Explications du problème

[Lien vers le sujet du projet (visité en fevrier 2021)](https://bthierry.pages.math.cnrs.fr/course-fem/projet/2020-2021/)

![le domaine](./2021-2021-flat.svg)

<img src="https://bthierry.pages.math.cnrs.fr/course-fem/_images/2020-2021-flat.svg">

Nous étudions le domaine $\Omega$ qui est un ouvert polygonal connexe. Les murs sont supposés parfaitement isolant, ce qui explique la condition de Neumann homogène que nous imposons. Nous remarquons que le bord de $\Omega$ est séparé en plusieurs parties : les radiateurs ($\Gamma_{\text{Rad}}$), les fenêtres ($\Gamma_{\text{Fen}}$) et les murs ($\Gamma_{\text{Mur}}$). Nous cherchons à calculer la température $u$ dans la pièce, qui vérifie le système suivant :
$$
(\mathrm{P}_\text{initial}) :
\left\{
\begin{array}{r c l l}
  -\Delta u & = & 0 & (\Omega) \\
  u & = & T_c & (\Gamma_{\text{Rad}})\\
  u & = & T_f & (\Gamma_{\text{Fen}})\\
  \partial_n u & = & 0 & (\Gamma_{\text{Mur}})
\end{array}
\right..
$$
Les paramètres sont les suivants :

- la longueur $L := 10$ ;
- la largeur $\ell := 10$ ;
- l'épaisseur des murs $d := 0.5$ ;
- la longueur d'une fenêtre est de 1 ;
- la longueur d'un radiateur est de 1 ;
- les températures $T_c := 25$ et $T_f := -10$ sont les températures respectivement des radiateurs et de dehors.

Nous souhaitons résoudre le problème $(\mathrm{P}_\text{initial})$ à l'aide de la méthode des éléments finis $\mathbb{P}^1$-Lagrange.

## Résolution du problème

### Formulation faible du problème

Considèrons l'espace de Hilbert $H^1 (\Omega)$ ainsi que $\gamma_{\Gamma_{\text{Rad}}}$ et $\gamma_{\Gamma_{\text{Fen}}}$ les applications trace de $H^1 (\Omega)$ sur $\mathrm{L}^2 (\Gamma_{\text{Rad}})$ et sur $\mathrm{L}^2 (\Gamma_{\text{Fen}})$ respectivement.
Pour résoudre $(\mathrm{P}_\text{initial})$, nous nous ramenons au cas de conditions de Dirichlet homogènes en introduisant un relèvement $u_\gamma \in H^1_{T_c, T_f} (\Omega)$ de $T_f$ et de $T_c$, avec
$$H^1_{T_c, T_f} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = T_c \wedge \gamma_{\Gamma_{\text{Fen}}} u = T_f \}.$$
De plus, nous définissons
$$H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = 0 \wedge \gamma_{\Gamma_{\text{Fen}}} u = 0 \}.$$
Le problème $(\mathrm{P}_\text{initial})$ revient alors à chercher $u_r := u - u_\gamma \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ satisfaisant :
$$
(\mathrm{P}_\text{relèvement}) :
\left\{
\begin{array}{r c l l}
  -\Delta u_r & = & 0 & (\Omega) \\
  u_r & = & 0 & (\Gamma_{\text{Rad}})\\
  u_r & = & 0 & (\Gamma_{\text{Fen}})\\
  \partial_n u_r & = & 0 & (\Gamma_{\text{Mur}})
\end{array}
\right..
$$
Multiplions la première équation de $(\mathrm{P}_\text{relèvement})$ par des fonctions tests $v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, intégrons sur $\Omega$ et appliquons le théorème de Green : **TODO laplacien de u_gamma**
$$
\begin{aligned}
  -\Delta u_r = 0
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ -\Delta u_r \cdot v = 0 \\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ -\int_{\Omega} \Delta u_r \cdot v  = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma} \partial_n u_r \cdot v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma_{\text{Rad}}} \partial_n u_r \cdot \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Fen}}} \partial_n u_r \cdot \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Mur}}} \underbrace{\partial_n u_r}_{=\ 0} \cdot v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u \cdot \nabla v - \int_{\Omega} \nabla u_\gamma \cdot \nabla v = 0 **TODO**
\end{aligned}
$$
Nous obtenons ainsi la formulation faible du problème :
$$
(\mathrm{P_{FF}}) :
\left\{
\begin{array}{l}
  \text{Trouver } u \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \text{ tel que }\\
  \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ a(u,v)=\ell(v)
\end{array}
\right.
$$
avec
$$
\begin{array}{r r c l}
  a :& H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \times H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) & \longrightarrow & \mathbb{R}\\
  & (u,v) & \longmapsto & \displaystyle \int_{\Omega}\nabla u\cdot\nabla v **TODO**\\
  \ell :& H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) & \longrightarrow & \mathbb{R}&&\\
  & v & \longmapsto & 0 **TODO**
\end{array}
$$

### Existence et unicité de la solution

Tentons d'appliquer le théorème de Lax-Milgram à $(\mathrm{P_{FF}})$.

1. $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ est un espace de Hilbert
2. $\ell(\cdot)$ est clairement linéaire (du fait de l'intégrale)
3. $a(\cdot,\cdot)$  est bilinéaire, pour la même raison
4. Continuité de $\ell(\cdot)$ : prenons une fonction $v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ : **TODO**
$$
\begin{aligned}
  | \ell(v) |
  & = \left| \int_{\Omega} fv \right|\\
  & \leqslant \| f \|_{\mathrm{L}^2(\Omega)} \| v \|_{\mathrm{L}^2(\Omega)} & \text{Cauchy-Schwarz}\\
  & \leqslant \| f \|_{\mathrm{L}^2(\Omega)} \| v \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} & \text{inégalité des normes} \\
\end{aligned}
$$
5. Continuité de $a(\cdot,\cdot)$ : prenons deux fonctions $u$ et $v$ de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ :
$$
\begin{aligned}
  | a(u, v) |
  & = \left| \int_{\Omega} \nabla u \cdot \nabla v \right|\\
  & \leqslant  \| \nabla u \|_{\mathrm{L}^2(\Omega)} \| \nabla v \|_{\mathrm{L}^2(\Omega)} & \text{inégalité triangulaire dans }  \mathrm{L}^2(\Omega)\\
  & \leqslant \| \nabla u \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} \| \nabla v \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} & \text{inégalité des normes} \\
\end{aligned}
$$
6. Coercivité de $a(\cdot, \cdot)$ : prenons une fonction $u \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ :
$$
\begin{aligned}
  a(u,u)
  & = \int_{\Omega} \nabla u \cdot \nabla u\\
  & = \int_{\Omega} \|\nabla u\|^2\\
  & = \| {u} \|^2_{\mathrm{L}^2(\Omega)}\\
  & \geqslant  C \| {u} \|^2_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} & \text{inégalité de Poincaré}
\end{aligned}
$$
Toutes les conditions sont réunies : le problème admet une unique solution d'après le théorème de Lax-Milgram.

### Maillage triangulaire (ou triangulation)

Nous découpons maintenant le domaine en triangles pour obtenir un maillage triangulaire (ou triangulation) conforme de $\Omega$. Ce découpage n'est pas fait manuellement, nous utilisons [GMSH](https://gmsh.info) ([un tutoriel](https://bthierry.pages.math.cnrs.fr/tutorial/gmsh) est fourni par Bertrand Thierry). Une telle triangulation sera notée $\mathcal{T}_h := \{K_p,\ p \in \llbracket 1, N_t \rrbracket\}$, l'indice $h$ faisant référence à la *finesse de maillage*, que l'on définit par le plus grand diamètre des triangles :
$$
h := \max_{K \in \mathcal{T}_h}(\mathrm{diam}\,(K)).
$$
Le diamètre d'un triangle est la distance maximale entre deux points du triangle. Nous notons de plus $\mathcal{S}_h$ l'ensemble des sommets de $\mathcal{T}_h$.

Nous rappelons qu'ici $\mathbb{P}^1$ est l'espace des polynômes réels de degré 1 sur $\omega \subset \mathbb{R}^2$ un ouvert, de dimension 3 :
$$\mathbb{P}^1 (\omega) := \{ p : \omega \to \mathbb{R} \ |\ \exists !a,b,c \in \mathbb{R}\ /\ \forall (x,y) \in \omega, p(x,y) = a + bx + cy \}.
$$
Nous pouvons maintenant introduire l'espace fonctionnel $\mathbb{P}^1$-Lagrange sur $\Omega$, souvent abrégé $\mathbb{P}^1$ et noté $V_h$. Il contient les fonctions continues sur $\overline{\Omega}$ et linéaires sur chaque triangle de $\mathcal{T}_h$ :
$$
V_h := \left\{ v_h \in \mathcal{C}^0 \left(\overline{\Omega}\right) \ |\ \forall K \in \mathcal{T}_h, v_h|_{K} \in \mathbb{P}^1(K) \right\}.
$$
En notant $N_S := \mathrm{card}\,(\mathcal{S}_h)$ le nombre de sommets du maillage $\mathcal{T}_h$, introduisons la famille des fonctions de forme $(\varphi_I)_{1\leqslant I \leqslant N_S}$ de $V_h$, qui sont nulles sur chaque sommet sauf un : le sommet $\mathrm{s}_I$. La famille $(\varphi_I)_{1\leqslant I \leqslant N_S}$ est une base de $V_h$, qui est alors de dimension $N_S$. Ainsi, nous pouvons écrire :
$$u_h = \sum_{I\ =\ 1}^{N_S} u_I\varphi_I,$$
avec $u_I = u_h(\mathrm{s}_I)$.

L'espace $V_h$ n'est pas tout à fait celui qui nous intéresse car nous allons montrer que nous voulons une *discrétisation* de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$. C'est pourquoi nous posons :
$$
V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} := \{ v_h \in V_h\ |\ v_h|_{\Gamma_{\text{Rad}}} = 0 \wedge \ v_h|_{\Gamma_{\text{Fen}}} = 0 \}.
$$
**TODO est-ce qu'on peut utiliser la même base ?**

### Méthode de Galerkin

Nous avons montré que le théorème de Lax-Milgram s'applique à la formulation variationnelle $(\mathrm{P_{FF}})$, et donc, que le problème $(\mathrm{P}_\text{initial})$ admet une unique solution. **TODO: vrai?**
Utilisons la méthode de Galerkin pour *approcher* l'espace $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ par un espace de Hilbert (pour le même produit scalaire) $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, de *dimension finie*. La formulation faible $(\mathrm{P_{FF}})$ est alors résolue dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$ uniquement. Ainsi, le problème se récrit :
$$
(\mathrm{P_{approché}}) :
\left\{
\begin{array}{l}
  \text{Trouver } u_h \in V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \text{ tel que}\\
  \forall v_h \in V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}},\ a(u_h, v_h) = \ell(v_h)
\end{array}
\right..
$$
Le problème *approché* $(\mathrm{P_{approché}})$ admet une unique solution. En effet, l'espace $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ est un sous-espace de Hilbert de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, nous pouvons donc appliquer le théorème de Lax-Milgram, dont les hypothèses sur $a(\cdot,\cdot)$ et $\ell(\cdot)$ sont toujours vérifiées dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$.

En élément finis $\mathbb{P}^1$, un relèvement de $T_c$ et de $T_f$ est la fonction $u_{\gamma}^h$ de $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$ telle que
$$
u_{\gamma}^h (\mathrm{s}_j) =
\left\{
\begin{array}{l l}
  T_c & \text{si }\mathrm{s}_j\in\Gamma_{\text{Rad}},\\
  T_f & \text{si }\mathrm{s}_j\in\Gamma_{\text{Fen}},\\
  0 & \text{sinon.}
\end{array}
\right.
$$
Cette fonction est un relèvement de l'interpolée de $T_c$ et de $T_f$ dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$.

### Formulation matricielle

Grâce à $(\mathrm{P_{approché}})$ et les propriétés des applications et de l'espace $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$, nous pouvons récrire le problème comme la résolution du système linéaire :
$$\mathrm{\bf A} \mathrm{\bf u}_h = \mathrm{\bf b}.$$
Les coefficients de la matrice $\mathrm{\bf A}$ et des vecteurs $\mathrm{\bf u}_h$ et $\mathrm{\bf b}$ sont donnés par :
$$
\begin{aligned}
    \mathrm{\bf A} & = (\mathrm{A}_{I,J})_{1 \leqslant I,J \leqslant N_S}, & \mathrm{A}_{I,J} &= a(\varphi_J,\varphi_J) = \int_{\Omega}\nabla \varphi_J\cdot\nabla\varphi_I **TODO**\\
    \mathrm{\bf u}_h & = (u_I)_{1 \leqslant I \leqslant N_S} & &\\
    \mathrm{\bf b} & = (\mathrm{b}_I)_{1 \leqslant I \leqslant N_S}, & \mathrm{b}_I & = \ell(\varphi_I) = 0 **TODO**
\end{aligned}
$$

Nous séparons les degrés de liberté en deux sous-ensembles (quitte à renuméroter) :

1. Ceux qui appartiennent à $\Omega$ ou à $\Gamma_{\text{Mur}}$ : nous les noterons avec un indice $I$ (pour Intérieur) ;
2. Ceux qui appartiennent à $\Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}$, ils seront notés avec un indice $D$. Le système linéaire devient :
$$
\mathrm{\bf A} \mathrm{\bf u}_h = \mathrm{\bf b} \iff \left(
\begin{array}{c c}
  \mathrm{\bf A}_{I,I}  & \mathrm{\bf A}_{I, D}\\
  \mathrm{\bf A}_{D, I} & \mathrm{\bf A}_{D,D}
\end{array}
\right) \left(
\begin{array}{c}
  \mathrm{\bf u}_I\\
  \mathrm{\bf u}_D
\end{array}
\right) =  \left(
\begin{array}{c}
  \mathrm{\bf b}_I\\
  \mathrm{\bf b}_D
\end{array}
\right)
$$
Nous notons $\mathrm{\bf u}_{\gamma}^h$ le vecteur de même taille que $\mathrm{\bf b}_D$ contenant les évaluation de $u_{\gamma}^h (\mathrm{s}_I)$ avec $\mathrm{s}_I \in \Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}$. Appliquer la condition de Dirichlet hétérogène se traduit par :
$$
\left(
\begin{array}{c c}
  \mathrm{\bf A}_{I,I} & \mathrm{\bf A}_{I,D}\\
  \mathrm{\bf 0} & \mathrm{\bf I}_{D,D}
\end{array}
\right)
\left(
\begin{array}{c}
  \mathrm{\bf u}_I\\
  \mathrm{\bf u}_D
\end{array}
\right)  =   \left(
\begin{array}{c}
  \mathrm{\bf b}_I\\
  \mathrm{\bf u}_{\gamma}^h
\end{array}
\right).
$$
La matrice obtenue est non symétrique, ce qui peut poser des problèmes (augmentation du coût de stockage mémoire, etc.). Une astuce simple consiste à récrire le système sous la forme suivante :
$$
\left(
\begin{array}{c c}
  \mathrm{\bf A}_{I,I} & \mathrm{\bf 0}\\
  \mathrm{\bf 0} & \mathrm{\bf I}_{D,D}
\end{array}
\right)
\left(
\begin{array}{c}
  \mathrm{\bf u}_I\\
  \mathrm{\bf u}_D
\end{array}
\right)  =   \left(
\begin{array}{c}
  \mathrm{\bf b}_I - \mathrm{\bf A}_{I,D} \mathrm{\bf u}_{\gamma}^h\\
  \mathrm{\bf u}_{\gamma}^h
\end{array}
\right).
$$
Nous pouvons aussi nous contenter de résoudre un système plus petit :
$$\mathrm{\bf A}_{I,I} \mathrm{\bf u}_I = \mathrm{\bf b}_I - \mathrm{\bf A}_{I,D} \mathrm{\bf u}_{\gamma}^h.$$

### Algorithme d'assemblage

Nous récrivons la matrice $\mathrm{\bf A}$ sous la forme suivante et calculons les *contributions élémentaires*, qui vont s'ajouter petit à petit dans la matrice $\mathrm{\bf A}$ :
$$
\begin{aligned}
  \mathrm{\bf A}
  & = \sum_{I\ =\ 1}^{N_s}\sum_{j\ =\ 0}^{N_s-1}a(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
  & = \sum_{p\ =\ 1}^{N_t}\sum_{I\ =\ 1}^{N_s}\sum_{J\ =\ 1}^{N_s} \underbrace{\int_{K_p}\nabla \varphi_J \cdot\nabla \varphi_I}_{\text{contribution élémentaire}} \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
\end{aligned}
$$
où $\mathrm{\bf e}_I$ est le vecteur de la base canonique de $\mathbb{R}^{N_s}$.

Nous devons maintenant travailler localement dans chaque triangle. Pour cela, nous avons besoin d'introduire une numérotation locale de chaque sommet une fonction $\mathrm{locToGlob}$ (*Local To Global*) permettant de basculer du local vers le global telle que, pour $p \in \llbracket 1,N_t \rrbracket$ et $i \in \llbracket 1,3 \rrbracket$ :
$$\mathrm{locToGlob}\,(p,i) = I \iff \mathrm{s}_i^p = \mathrm{s}_I.$$
Nous distinguerons la numérotation globale par des lettres capitales ($I$, $J$) et la numérotation locale par des minuscules ($i$, $j$). Nous introduisons aussi les fonctions de forme locales :
$$ \varphi_i^p = \varphi_{\mathrm{locToGlob}\,(p,i)}|_{K_p}.$$
Utilisons ces nouvelles notations en ramenant la somme sur tous les sommets du maillage à uniquement les sommets du triangle considéré :
$$
\mathrm{\bf A} = \sum_{p=1}^{N_t}\sum_{i=1}^{3}\sum_{j=1}^{3} \int_{K_p}\nabla \varphi_j^p \cdot\nabla \varphi_i^p\ \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}^\top\mathrm{\bf e}_{\mathrm{locToGlob}\,(p,j)}.
$$

Le pseudo-code de l'algorithme d'assemblage :

```
A = 0
b = 0
For p = 1:N_t
  For i = 1:3
    I = locToGlob(p,i)
    For j = 1:3
      J = locToGlob(p,j)
      A(I,J) += ∫_{K_p}(∇ϕ_j^p·∇ϕ_i^p)
    EndFor
    b(I) += l_p(ϕ_i^p)
  EndFor
EndFor
```

<!-- #### Calcul des matrices élémentaires **TODO**

Dans le triangle de référence $\widehat{K}$, la matrice de rigidité élémentaire $\mathrm{\bf \widehat{D}}$ a pour expression
$$
\mathrm{\bf \widehat{D}} = \frac{1}{2}
\left(
\begin{array}{l l c}
  2 & -1 & -1 \\
  -1 & 1 & 0 \\
  -1 & 0 & 1
\end{array}
\right)
$$
Les coefficients de la matrice de rigidité élémentaire $\mathrm{\bf D}^e_p = ((\mathrm{D}^e{p})_{i,j})_{1\leq i,j\leq 3}$ sont obtenus pas la relation suivante :
$$
\begin{aligned}
  (\mathrm{D}^e_p)_{i,j} &= \int_{K_p}\nabla \varphi_j^p \cdot \nabla\varphi_i^p,\\
  &= | K_p | (\nabla\widehat{\varphi}_j)^T \left( \left(\left(J_p^\top \right)^{-1} \right)^\top \left( J_p^\top \right)^{-1} \right) \nabla\widehat{\varphi}_i. **TODO**
\end{aligned}
$$ -->

Les coefficients sont obtenus pas la relation suivante :
$$
\int_{K_p}\nabla \varphi_j^p \cdot \nabla\varphi_i^p = | K_p | (\nabla\widehat{\varphi}_j)^\top \left( \left(\left(J_p^\top \right)^{-1} \right)^\top \left( J_p^\top \right)^{-1} \right) \nabla\widehat{\varphi}_i.
$$
 **TODO**
avec $J_p$ la matrice jacobienne de la tranformation des coordonnées de $\widehat{K}$ en celles de $K_p$ et où les gradients des fonctions de forme sur $\widehat{K}$ $\widehat{\varphi}_i$ sont donnés par :
$$
\nabla\widehat{\varphi}_0 =
\begin{pmatrix}
  -1\\
  -1
\end{pmatrix}
,
\quad
\nabla\widehat{\varphi}_1 =
\begin{pmatrix}
  1\\
  0
\end{pmatrix},
\quad
\nabla\widehat{\varphi}_2 =
\begin{pmatrix}
  0\\
  1
\end{pmatrix}.
$$

## Implantation

Utilisation des format *COO* et *CSR* pour stocker les matrices et faire les calculs.

Informatiquement, nous devons donc rendre les lignes et colonnes associées aux degrés de liberté de Dirichlet, nulles, sauf sur la diagonale avec la valeur 1. Cette opération peut être effectuée après l'assemblage de la matrice ou lors de l'algorithme directement.