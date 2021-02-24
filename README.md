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
  - [ ] matrice de rigidité
  - [ ] assemblage des matrices élémentaires
- 3
  - [ ] faire `main.py`
  - [ ] affichage graphique avec gradient de couleur
- 4
  - [ ] commenter le code
  - [ ] faire le readme : application trace, définir le bon Vh, vérifier la formulation, formulation matricielle

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
(P_{initial}) :
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

Nous souhaitons résoudre le problème $(P_{initial})$ à l'aide de la méthode des éléments finis $\mathbb{P}^1$-Lagrange.

## Résolution du problème

### Formulation faible du problème

Nous nous ramenons au cas de conditions de Dirichlet homogènes en introduisant un *relèvement* $u_\gamma \in H^1_{T_c, T_f} (\Omega)$ de $T_f$ et de $T_c$, avec

$$H^1_{T_c, T_f} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = T_c,\ \gamma_{\Gamma_{\text{Fen}}} u = T_f \}.$$

De plus, définissons

$$H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = 0,\ \gamma_{\Gamma_{\text{Fen}}} u = 0 \}.$$

Le problème $(P_{initial})$ revient alors à chercher $u_r := u - u_\gamma \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ satisfaisant :

$$
(P_{relèvement}) :
\left\{
\begin{array}{r c l l}
  -\Delta u_r & = & 0 & (\Omega) \\
  u_r & = & 0 & (\Gamma_{\text{Rad}})\\
  u_r & = & 0 & (\Gamma_{\text{Fen}})\\
  \partial_n u_r & = & 0 & (\Gamma_{\text{Mur}})
\end{array}
\right..
$$

Multiplions la première équation de $(P_{relèvement})$ par des fonctions tests $v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, intégrons sur $\Omega$ et appliquons le théorème de Green :

$$
\begin{aligned}
  -\Delta u_r = 0
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ (-\Delta u_r) v = 0 \\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ -\int_{\Omega}(\Delta u_r) v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma} (\partial_n u_r) v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma_{\text{Rad}}} (\partial_n u_r) \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Fen}}} (\partial_n u_r) \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Mur}}} (\underbrace{\partial_n u_r}_{=\ 0}) v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u_r \cdot \nabla v = 0\\
  & \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ \int_{\Omega} \nabla u \cdot \nabla v - \int_{\Omega} \nabla u_\gamma \cdot \nabla v = 0 TODO
\end{aligned}
$$

Nous obtenons ainsi la formulation faible du problème :

$$
(P_{FF}) :
\left\{
\begin{array}{l}
  \text{Trouver } u \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \text{ tel que } TODO \\
  \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ a(u,v)=\ell(v)
\end{array}
\right.
$$

avec $a(\cdot,\cdot) : H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \times H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)  \to \mathbb{R}$ et $\ell(\cdot) : H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \to \mathbb{R}$ définies par

$$
\left\{
\begin{aligned}
  a(u,v) & := \int_{\Omega} \nabla u \cdot \nabla v - \int_{\Omega} \nabla u_\gamma \cdot \nabla v TODO \\
  \ell(v) & := 0
\end{aligned}
\right..
$$

___

### Existence et unicité de la solution TODO

Tentons d'appliquer le théorème de Lax-Milgram à cette formulation faible

1. $\Ho$ est un espace de Hilbert
2. $\ell(\cdot)$ est clairement linéaire (du fait de l'intégrale)
3. $a(\cdot,\cdot)$  est bilinéaire, pour la même raison
4. Continuité de $\ell(\cdot)$ : prenons une fonction $v\in\Ho$ :

$$
\begin{aligned}
  \abs{\ell(v)}
  & = \underbrace{\abs{\int_{\Omega} fv}}_{\PSL{f}{v}}\\
  & \leq  \normL{f}\normL{v} & \text{Cauchy-Schwarz}\\
  & \leq   \underbrace{\normL{f}}_{\text{Constant}}\normH{v} & \text{inégalité des normes} \\
\end{aligned}
$$

5. Continuité de $a(\cdot,\cdot)$ : prenons deux fonctions $u$ et $v$ de $\Ho$ :

$$
\begin{aligned}
  \abs{a(u,v)}
  &= \abs{\int_{\Omega} \nabla u \cdot \nabla v + c\int_{\Omega} u v}\\ 
  & \leq  \underbrace{\abs{\int_{\Omega} \nabla u \cdot \nabla v}}_{\PSLd{\nabla u}{\nabla v}} + \abs{c}\underbrace{\abs{\int_{\Omega} u v}}_{\PSL{u}{v}} & \text{inégalité classique}\\
  & \leq  \normLd{\nabla u}\normLd{\nabla v} + \abs{c} \normL{u}\normL{v} & \text{inégalité triangulaire dans}  \Lo\\
  & \leq   \normH{u}\normH{v}+ \abs{c} \normH{u}\normH{v} & \text{inégalité des normes} \\
  & \leq   (1+c)\normH{u}\normH{v}
\end{aligned}
$$

6. Coercivité de $a(\cdot, \cdot)$ : prenons une fonction $u\in\Ho$ :

$$
\begin{aligned}
  a(u,u)
  & = \int_{\Omega} \nabla u \cdot \nabla u + c\int_{\Omega} u u = \int_{\Omega} \|\nabla u\|^2 + c\int_{\Omega} |u|^2\\ 
  & \geq \min(1,c)\left(\int_{\Omega} \|\nabla u\|^2 + \int_{\Omega} |u|^2\right)\\ 
  & \geq \min(1,c)\normH{u}^2 
\end{aligned}
$$

Toutes les conditions sont réunies : le problème admet une unique solution d'après le théorème de Lax-Milgram.

___

### Méthode de Galerkin

Nous considérons la formulation variationnelle $(P_{FF})$. Nous avons montré que le théorème de Lax-Milgram s'applique, et donc, que le problème $(P_{initial})$ admet une unique solution. TODO
Utilisons la méthode de Galerkin pour *approcher* l'espace $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ par un espace de Hilbert, pour le même produit scalaire, $V_h \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, de *dimension finie*. La formulation faible $(P_{FF})$ est alors résolue dans $V_h$ uniquement, avec pour solution $u_h$ :

$$
(P_{approché})
\left\{
\begin{array}{l}
  \text{Trouver } u_h \in V_h \text{ tel que}\\
  \forall v_h \in V_h,\ a(u_h, v_h) = \ell(v_h).
\end{array}
\right.
$$

Le problème *approché* $(P_{approché})$ admet une unique solution. En effet, L'espace $V_h \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ est un sous-espace de Hilbert de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, nous pouvons donc appliquer le théorème de Lax-Milgram, dont les hypothèses sur $a(\cdot,\cdot)$ et $\ell(\cdot)$ sont toujours vérifiées dans $V_h$.

#### Maillage triangulaire (ou triangulation)

Nous découpons maintenant le domaine en triangles pour obtenir un maillage triangulaire (ou triangulation) conforme de $\Omega$. Une telle triangulation sera notée $\mathcal{T}_h := \{K_p,\ p \in \llbracket 1, N_t \rrbracket\}$, l'indice $h$ faisant référence à la *finesse de maillage*, que l'on définit par le plus grand diamètre des triangles :

$$
h := \max_{K \in \mathcal{T}_h}(\mathrm{diam}\,(K)).
$$

Le diamètre d'un triangle est la distance maximale entre deux points du triangle. Nous notons de plus $\mathcal{S}_h$ et $\mathcal{A}_h$ les ensembles respectivement des sommets et des arêtes de $\mathcal{T}_h$.

Ce découpage n'est pas fait manuellement, nous utilisons [GMSH](https://gmsh.info) ([un tutoriel](https://bthierry.pages.math.cnrs.fr/tutorial/gmsh) est fourni par Bertrand Thierry).

Nous rappelons qu'ici $\mathbb{P}^1$ est l'espace des polynômes réels de degré 1 sur $\omega \subset \mathbb{R}^2$ un ouvert, de dimension 3 :

$$\mathbb{P}^1 (\omega) := \{ p : \omega \to \mathbb{R} \ |\ \exists !a,b,c \in \mathbb{R}\ /\ \forall (x,y) \in \omega, p(x,y) = a + bx + cy \}.
$$

Nous pouvons maintenant introduire l'espace fonctionnel $\mathbb{P}^1$-Lagrange sur $\Omega$, souvent abrégé $\mathbb{P}^1$ et noté $V_h$. Il contient les fonctions continues sur $\overline{\Omega}$ et linéaires sur chaque triangle de $\mathcal{T}_h$ :

$$
V_h := \{ v_h \in \mathcal{C}^0(\overline{\Omega})\ |\ \forall K \in \mathcal{T}_h, v_h|_{K} \in \mathbb{P}^1(K) \}.
$$

En notant $N_S = \mathrm{card}\,(\mathcal{S}_h)$ le nombre de sommets du maillage $\mathcal{T}_h$, introduisons la famille des fonctions de forme $(\varphi_I)_{1\leqslant I \leqslant N_S}$ de $V_h$, qui sont nulles sur chaque sommet sauf un : le sommet $\mathrm{s}_I$. La famille $(\varphi_I)_{1\leqslant I \leqslant N_S}$ est une base de $V_h$, qui est alors de dimension $N_S$.

#### Formulation matricielle TODO

TODO

La matrice $A$ peut être décomposée en deux matrices :

$$A := D + c M$$

avec
$M$ la matrice de masse (ou de volume), de coefficient

  $$M_{I,J} = \int_{\Omega} \varphi_J\varphi_I.,$$

$D$ la matrice de rigidité, de coefficient

  $$D_{I,J}=  \int_{\Omega}\nabla\varphi_J \cdot \nabla\varphi_I.$$

Dans la littérature, cette matrice est souvent notée $K$, mais nous l'appelons $D$ pour éviter toute confusion avec les triangles, nommés $K$ également.

Ici, la matrice $A$ se récrit simplement : $A = D$

#### Algorithme d'assemblage

Nous parcourrons les triangles du maillage et calculons les *contributions élémentaires*, qui vont s'ajouter petit à petit dans la matrice $A$.

Reprenons l'expression du coefficient $A_{I,J}$:

$$
A_{I,J} = \int_{\Omega}\nabla \varphi_J \cdot\nabla \varphi_I = \sum_{p\ =\ 1}^{N_t} \underbrace{\int_{K_p}\nabla \varphi_J \cdot\nabla \varphi_I}_{\text{contribution élémentaire}} TODO.
$$

Introduisons $a_p(\cdot,\cdot)$ la famille de forme bilinéaire suivante, pour $p \in \llbracket 1,N_t \rrbracket$ : 

$$
a_p(\varphi_J,\varphi_I) = \int_{K_p}\nabla \varphi_J \cdot \nabla \varphi_I.
$$
Ensuite, nous réécrivons la matrice $A$ sous la forme suivante

$$
A = \sum_{I\ =\ 1}^{N_s}\sum_{j=0}^{N_s-1}a(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J,
$$
où $\mathrm{\bf e}_I$ est le vecteur de la base canonique de $\mathbb{R}^{N_s}$.  Nous avons alors

$$ 
\begin{aligned}
  A
  & = \sum_{I\ =\ 1}^{N_s}\sum_{J\ =\ 1}^{N_s}a(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
  & = \sum_{I\ =\ 1}^{N_s}\sum_{J\ =\ 1}^{N_s}\sum_{p\ =\ 1}^{N_t}a_{p}(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
  & = \sum_{p\ =\ 1}^{N_t}\sum_{I\ =\ 1}^{N_s}\sum_{J\ =\ 1}^{N_s}a_{p}(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
\end{aligned}
$$
Nous remarquons maintenant que $a_{p}(\varphi_J,\varphi_I)$ est nul dès lors que $\mathrm{s}_I$ ou $\mathrm{s}_J$ ne sont pas des sommets de $K_p$ (car $\varphi_I\varphi_J = 0$ sur $K_p$). Finalement, la somme sur tous les sommets du maillage se réduit à une somme sur les 3 sommets du triangle $K_p$ considéré. 

Nous comprenons que nous devons maintenant travailler localement dans chaque triangle. Pour cela, nous avons besoin d'introduire une numérotation locale de chaque sommet une fonction $\mathrm{locToGlob}$ (*Local To Global*) permettant de basculer du local vers le global une fonction telle que, pour $p \in \llbracket 1,N_t \rrbracket$ et $i \in \llbracket 1,3 \rrbracket$ :

$$\mathrm{locToGlob}\,(p,i) = I \iff \mathrm{s}_i^p = \mathrm{s}_I.$$

Ainsi, pour un triangle  $K_p$, ses sommets sont numérotés $[\mathrm{s}_{1}^{p},\mathrm{s}_{2}^{p},\mathrm{s}_{3}^{p}]$ en numérotation locale ou $[\mathrm{s}_{\mathrm{locToGlob}\,(p,1)},\mathrm{s}_{\mathrm{locToGlob}\,(p,2)},\mathrm{s}_{\mathrm{locToGlob}\,(p,3)}]$ en numérotation globale, comme le montre la figure :numref:${number} <fig-loc2glob>$. Nous distinguerons la numérotation globale par des lettres capitales ($I$, $J$) et la numérotation locale par des minuscules ($i$, $j$). Nous introduisons aussi les fonctions de forme locales :

$$ \varphi_i^p = \varphi_{\mathrm{locToGlob}\,(p,i)}|_{K_p}.
$$
Utilisons ces nouvelles notations dans l'équation :eq:$eq-assemble_tmp$, en ramenant la somme sur les sommets à uniquement les sommets du triangle considéré :

$$ A = \sum_{p=1}^{N_t}\sum_{i=1}^{3}\sum_{j=1}^{3}a_{p}(\varphi_j^p,\varphi_i^p) \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}^\top\mathrm{\bf e}_{\mathrm{locToGlob}\,(p,j)}
$$

Le  pseudo-code de l'algorithme d'assemblage
```
A = 0
B = 0
For p = 1:N_t
  For i = 1:3
    I = L2G(p,i)
    For j = 1:3
      J = L2G(p,j)
      A(I,J) += a_p(ϕ_j^p,ϕ_i^p)
    EndFor
    B(I) += l_p(ϕ_i^p)
  EndFor
EndFor
```

Les *contributions élémentaires*, c'est à dire les quantités :math:`a_p(\mphi_j^p,\mphi_i^p)` et :math:`\ell_{p}(\mphi_i^p)`, peuvent elles aussi être décomposées en deux parties. Pour rappel, les sommets d'un triangle :math:`\tri_p` seront notés :math:`[\vertice_{0}^{p}, \vertice_{1}^{p},\vertice_{2}^{p}]` et ordonnés dans le sens trigonométrique. Nous noterons :math:`\vertice_i^p=(x_i^p, y_i^p)` un sommet de :math:`\tri_p` et :math:`\mphi_i^p` la fonction de forme locale associée. Nous notons :math:`\De{p}` la matrice de rigidité élémentaire du triangle :math:`\tri_p`, de coefficient $(\De{p})_{i,j} &=\int_{\tri_p}\nabla\mphi_j^p\cdot\nabla\mphi_i^p.$


#### Calcul des matrices élémentaires
