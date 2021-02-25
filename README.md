# Projet de maillage et éléments finis

par Jérôme Bonacchi et Homer Durand à Polytech Sorbonne en spécialité mathématiques appliquées et informatique

- [Projet de maillage et éléments finis](#projet-de-maillage-et-éléments-finis)
  - [TODO](#todo)
  - [Exécution](#exécution)
  - [Explications du problème](#explications-du-problème)
  - [Résolution du problème](#résolution-du-problème)
    - [Formulation faible du problème](#formulation-faible-du-problème)
    - [Existence et unicité de la solution](#existence-et-unicité-de-la-solution)
    - [Maillage triangulaire](#maillage-triangulaire)
    - [Méthode de Galerkin](#méthode-de-galerkin)
    - [Formulation matricielle](#formulation-matricielle)
    - [Algorithme d'assemblage](#algorithme-dassemblage)
    - [Calcul des *contributions élémentaires*](#calcul-des-contributions-élémentaires)
  - [Implantation](#implantation)

## TODO

- 1a
  - [x] gmshToMesh
  - [x] getElements, getPoints
- 1b
  - [x] modélisation en gmsh
- 1c
  - [x] récrire le problème/la reformulation au propre
- 2a
  - [x] tester la modélisation gmsh : segments ne contient que les aretes des groupes physiques de dimension 1, si on veut les aretes des triangles ils faut réfléchir car les aretes des triangles ne sont pas des éléments et des entités en elles-mêmes. Bref, on s'en fout, non ?
- 2b
  - [x] matrice de rigidité
- 2c
  - [ ] ~~calculer la quadrature du membre de droite~~
- 2d
  - [x] faire l'assemblage des matrices élémentaires
- 2e
  - [x] faire locToGlob
- 2f
  - [x] condition de Dirichlet
- 2g
  - [ ] implémenter le calcul du membre de droite
- 3a
  - [ ] faire `main.py`
  - [x] vérifier la matrice de rigidité avec DU=0
  - [ ] vérifier locToGlob
  - [ ] vérifier condition de Dirichlet
  - [ ] affichage graphique avec gradient de couleur
- 4a
  - [ ] commenter le code
- 4b
  - [ ] formater le code
- 4c
  - [ ] faire le readme : s'occuper des TODO, expliquer l'implantation

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
$$H^1_{T_c, T_f} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = T_c \land \gamma_{\Gamma_{\text{Fen}}} u = T_f \}.$$
De plus, nous définissons
$$H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) := \{ u\in H^1 (\Omega)\ |\ \gamma_{\Gamma_{\text{Rad}}} u = 0 \land \gamma_{\Gamma_{\text{Fen}}} u = 0 \}.$$
Le problème $(\mathrm{P}_\text{initial})$ revient alors à chercher $u_r := u - u_\gamma \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ satisfaisant :
$$
(\mathrm{P}_\text{relèvement}) :
\left\{
\begin{array}{r c l l}
  -\Delta u_r & = & \Delta u_\gamma & (\Omega) \\
  u_r & = & 0 & (\Gamma_{\text{Rad}})\\
  u_r & = & 0 & (\Gamma_{\text{Fen}})\\
  \partial_n u_r & = & - \partial_n u_\gamma & (\Gamma_{\text{Mur}})
\end{array}
\right..
$$
Multiplions la première équation de $(\mathrm{P}_\text{relèvement})$ par des fonctions tests $v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, intégrons sur $\Omega$ et appliquons le théorème de Green :
$$
\begin{aligned}
  & -\Delta u_r = \Delta u_\gamma\\
  \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),& -\Delta u_r \cdot v = \Delta u_\gamma \cdot v \\
  \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),& -\int_{\Omega} \Delta u_r \cdot v  = \int_\Omega \Delta u_\gamma \cdot v\\
  \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),& \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma} \partial_n u_r \cdot v = - \int_{\Omega} \nabla u_\gamma \cdot \nabla v + \int_{\Gamma} \partial_n u_\gamma \cdot v\\
  \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),& \int_{\Omega} \nabla u_r \cdot \nabla v - \int_{\Gamma_{\text{Rad}}} \partial_n u_r \cdot \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Fen}}} \partial_n u_r \cdot \underbrace{v}_{=\ 0} - \int_{\Gamma_{\text{Mur}}} \underbrace{\partial_n u_r}_{=\ - \partial_n u_\gamma} \cdot v = - \int_{\Omega} \nabla u_\gamma \cdot \nabla v + \int_{\Gamma_{\text{Rad}}} \partial_n u_\gamma \cdot \underbrace{v}_{=\ 0} + \int_{\Gamma_{\text{Fen}}} \partial_n u_\gamma \cdot \underbrace{v}_{=\ 0} + \int_{\Gamma_{\text{Mur}}} \partial_n u_\gamma \cdot v\\
  \implies \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),& \int_{\Omega} \nabla u_r \cdot \nabla v = - \int_{\Omega} \nabla u_\gamma \cdot \nabla v
\end{aligned}
$$
Nous obtenons ainsi la formulation faible du problème :
$$
(\mathrm{P_{FF}}) :
\left\{
\begin{array}{l}
  \text{Trouver } u_r \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \text{ tel que }\\
  \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ a(u_r,v)=\ell(v)
\end{array}
\right.
$$
avec
$$
\begin{array}{r r c l}
  a :& H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \times H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) & \longrightarrow & \mathbb{R}\\
  & (u,v) & \longmapsto & \displaystyle \int_{\Omega}\nabla u\cdot\nabla v\\
  \ell :& H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) & \longrightarrow & \mathbb{R}&\\
  & v & \longmapsto & \displaystyle - \int_{\Omega} \nabla u_\gamma \cdot \nabla v
\end{array}
$$

### Existence et unicité de la solution

Tentons d'appliquer le théorème de Lax-Milgram à $(\mathrm{P_{FF}})$.

1. $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ est un espace de Hilbert
2. $\ell(\cdot)$ est clairement linéaire (du fait de l'intégrale)
3. $a(\cdot,\cdot)$  est bilinéaire, pour la même raison
4. Continuité de $\ell(\cdot)$ :
$$
\begin{aligned}
  \forall v \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega),\ | \ell(v) |
  & = \left| -\int_{\Omega} \nabla u_\gamma \cdot \nabla v \right|\\
  & \leqslant \| \nabla u_\gamma \|_{\mathrm{L}^2(\Omega)} \| \nabla v \|_{\mathrm{L}^2(\Omega)} & \text{inégalité de Cauchy-Schwarz dans } \mathrm{L}^2(\Omega)\\
  & \leqslant \| \nabla u_\gamma \|_{\mathrm{L}^2(\Omega)} \| v \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} & \text{inégalité des normes} \\
\end{aligned}
$$
5. Continuité de $a(\cdot,\cdot)$ : $\forall (u, v) \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega) \times H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$,
$$
\begin{aligned}
  | a(u, v) |
  & = \left| \int_{\Omega} \nabla u \cdot \nabla v \right|\\
  & \leqslant  \| \nabla u \|_{\mathrm{L}^2(\Omega)} \| \nabla v \|_{\mathrm{L}^2(\Omega)} & \text{inégalité de Cauchy-Schwarz dans } \mathrm{L}^2(\Omega)\\
  & \leqslant \| \nabla u \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} \| \nabla v \|_{H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)} & \text{inégalité des normes} \\
\end{aligned}
$$
6. Coercivité de $a(\cdot, \cdot)$ : $\forall u \in H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$,
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

### Maillage triangulaire

Nous découpons maintenant le domaine en triangles pour obtenir un maillage triangulaire conforme de $\Omega$. Ce découpage n'est pas fait manuellement, nous utilisons [GMSH](https://gmsh.info). Une telle triangulation sera notée $\mathcal{T}_h := \{K_p,\ p \in \llbracket 1, N_t \rrbracket\}$, l'indice $h$ faisant référence à la *finesse du maillage*, que l'on définit par le plus grand diamètre des triangles :
$$
h := \max_{K \in \mathcal{T}_h}(\mathrm{diam}\,(K)).
$$
Le diamètre d'un triangle est la distance maximale entre deux points du triangle. Nous notons de plus $\mathcal{A}_h$ et $\mathcal{S}_h$ les ensembles respectifs des arêtes et des sommets de $\mathcal{T}_h$.

Nous rappelons qu'ici $\mathbb{P}^1$ est l'espace des polynômes réels de degré 1 sur $\omega \subset \mathbb{R}^2$ un ouvert :
$$\mathbb{P}^1 (\omega) := \{ p : \omega \to \mathbb{R} \ |\ \exists !a,b,c \in \mathbb{R}\ /\ \forall (x,y) \in \omega, p(x,y) = a + bx + cy \}.
$$
Nous pouvons maintenant introduire l'espace fonctionnel $\mathbb{P}^1$-Lagrange sur $\Omega$, souvent abrégé $\mathbb{P}^1$ et noté $V_h$. Il contient les fonctions continues sur $\overline{\Omega}$ et linéaires sur chaque triangle de $\mathcal{T}_h$ :
$$
V_h := \left\{ v_h \in \mathcal{C}^0 \left(\overline{\Omega}\right) \ |\ \forall K \in \mathcal{T}_h, v_h|_{K} \in \mathbb{P}^1(K) \right\}.
$$
En notant $N_S := \mathrm{card}\,(\mathcal{S}_h)$ le nombre de sommets du maillage $\mathcal{T}_h$, introduisons la famille des fonctions de forme $(\varphi_I)_{1\leqslant I \leqslant N_S}$ de $V_h$, qui sont nulles sur chaque sommet sauf un : le sommet $\mathrm{s}_I$. La famille $(\varphi_I)_{1\leqslant I \leqslant N_S}$ est une base de $V_h$, qui est alors de dimension $N_S$. Ainsi, nous pouvons écrire :
$$\forall u_h \in V_h, \quad u_h = \sum_{I\ =\ 1}^{N_S} u_I\varphi_I, \quad \text{avec} \quad u_I = u_h(\mathrm{s}_I).$$
L'espace $V_h$ n'est pas tout à fait celui qui nous intéresse car nous allons montrer que nous voulons une *discrétisation* de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$. C'est pourquoi nous posons :
$$
V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} := \{ v_h \in V_h\ |\ v_h|_{\Gamma_{\text{Rad}}} = 0 \land \ v_h|_{\Gamma_{\text{Fen}}} = 0 \}.
$$

### Méthode de Galerkin

Nous avons montré que le théorème de Lax-Milgram s'applique à la formulation variationnelle $(\mathrm{P_{FF}})$, et donc, que le problème $(\mathrm{P}_\text{initial})$ admet une unique solution. **TODO: vrai?**
Utilisons la méthode de Galerkin pour *approcher* l'espace $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ par un espace de Hilbert (pour le même produit scalaire) $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, de *dimension finie*. La formulation faible $(\mathrm{P_{FF}})$ est alors résolue dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$. Ainsi, le problème se récrit :
$$
(\mathrm{P_{approché}}) :
\left\{
\begin{array}{l}
  \text{Trouver } u_r^h \in V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \text{ tel que}\\
  \forall v_h \in V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}},\ a(u_r^h, v_h) = \ell(v_h)
\end{array}
\right..
$$
Le problème *approché* $(\mathrm{P_{approché}})$ admet une unique solution. En effet, l'espace $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} \subset H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$ est un sous-espace de Hilbert de $H^1_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}} (\Omega)$, nous pouvons donc appliquer le théorème de Lax-Milgram, dont les hypothèses sur $a(\cdot,\cdot)$ et $\ell(\cdot)$ sont toujours vérifiées dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$.

En élément finis $\mathbb{P}^1$, un relèvement de $T_c$ et de $T_f$ est la fonction $u_{\gamma}^h$ de $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$ telle que
$$
u_{\gamma}^h (\mathrm{s}_J) =
\left\{
\begin{array}{l l}
  T_c & \text{si }\mathrm{s}_J\in\Gamma_{\text{Rad}},\\
  T_f & \text{si }\mathrm{s}_J\in\Gamma_{\text{Fen}},\\
  0 & \text{sinon.}
\end{array}
\right.
$$
Cette fonction est un relèvement de l'interpolée de $T_c$ et de $T_f$ dans $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$.

### Formulation matricielle

Grâce à $(\mathrm{P_{approché}})$, les propriétés des applications $a(\cdot,\cdot)$ et $\ell(\cdot)$ et de l'espace $V^h_{\Gamma_{\text{Rad}}, \Gamma_{\text{Fen}}}$, nous pouvons récrire ce problème comme la résolution du système linéaire :
$$
(\mathrm{P_{matriciel}}) : \mathrm{\bf A} \mathrm{\bf u}_r^h = \mathrm{\bf b}.
$$
Les coefficients de la matrice $\mathrm{\bf A}$ et des vecteurs $\mathrm{\bf u}_r^h$ et $\mathrm{\bf b}$ sont donnés par :
$$
\begin{aligned}
    \mathrm{\bf A} & = (\mathrm{A}_{I,J})_{1 \leqslant I,J \leqslant N_S}, & \mathrm{A}_{I,J} &= a(\varphi_J,\varphi_J) = \int_{\Omega}\nabla \varphi_J\cdot\nabla\varphi_I\\
    \mathrm{\bf u}_r^h & = ({u_r^h}_I)_{1 \leqslant I \leqslant N_S} & &\\
    \mathrm{\bf b} & = (\mathrm{b}_I)_{1 \leqslant I \leqslant N_S}, & \mathrm{b}_I & = \ell(\varphi_I) = - \int_{\Omega} \nabla u_\gamma^h \cdot \nabla \varphi_I
\end{aligned}
$$
Nous séparons les degrés de liberté en deux sous-ensembles (quitte à renuméroter) :

1. Ceux qui appartiennent à $\Omega \cup \Gamma_{\text{Mur}}$ : nous les notons avec un indice $\mathcal{I}$ ;
2. Ceux qui appartiennent à $\Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}$, ils sont notés avec un indice $\mathcal{D}$.

Le système $(\mathrm{P_{matriciel}})$ devient :
$$
\mathrm{\bf A} \mathrm{\bf u}_r^h = \mathrm{\bf b} \iff \left(
\begin{array}{c c}
  \mathrm{\bf A}_{\mathcal{I},\mathcal{I}}  & \mathrm{\bf A}_{\mathcal{I}, \mathcal{D}}\\
  \mathrm{\bf A}_{\mathcal{D}, \mathcal{I}} & \mathrm{\bf A}_{\mathcal{D},\mathcal{D}}
\end{array}
\right) \left(
\begin{array}{c}
  {\mathrm{\bf u}_r^h}_\mathcal{I}\\
  {\mathrm{\bf u}_r^h}_\mathcal{D}
\end{array}
\right) =  \left(
\begin{array}{c}
  \mathrm{\bf b}_\mathcal{I}\\
  \mathrm{\bf b}_\mathcal{D}
\end{array}
\right)
$$
Appliquer la condition de Dirichlet homogène se traduit par :
$$
\left(
\begin{array}{c c}
  \mathrm{\bf A}_{\mathcal{I},\mathcal{I}} & \mathrm{\bf A}_{\mathcal{I},\mathcal{D}}\\
  \mathrm{\bf 0} & \mathrm{\bf \mathcal{I}}_{\mathcal{D},\mathcal{D}}
\end{array}
\right)
\left(
\begin{array}{c}
  {\mathrm{\bf u}_r^h}_\mathcal{I}\\
  {\mathrm{\bf u}_r^h}_\mathcal{D}
\end{array}
\right)  =   \left(
\begin{array}{c}
  \mathrm{\bf b}_\mathcal{I}\\
  \mathrm{\bf 0}
\end{array}
\right).
$$
Le système $(\mathrm{P_{matriciel}})$ se simplifie alors en :
$$
(\mathrm{P'_{matriciel}}) :
\mathrm{\bf A}_{\mathcal{I},\mathcal{I}} {\mathrm{\bf u}_r^h}_\mathcal{I} = \mathrm{\bf b}_\mathcal{I}.
$$

### Algorithme d'assemblage

Nous récrivons la matrice $\mathrm{\bf A}$ sous la forme suivante et calculons les *contributions élémentaires*, qui vont s'ajouter petit à petit dans la matrice $\mathrm{\bf A}$ :
$$
\begin{aligned}
  \mathrm{\bf A}
  & = \sum_{I\ =\ 1}^{N_s}\sum_{j\ =\ 0}^{N_s-1}a(\varphi_J,\varphi_I) \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
  & = \sum_{p\ =\ 1}^{N_t}\sum_{I\ =\ 1}^{N_s}\sum_{J\ =\ 1}^{N_s} \underbrace{\int_{K_p}\nabla \varphi_J \cdot\nabla \varphi_I}_{\text{contribution élémentaire}} \mathrm{\bf e}_I^\top\mathrm{\bf e}_J\\
\end{aligned}
$$
où $\mathrm{\bf e}_I$ est le vecteur de la base canonique de $\mathbb{R}^{N_s}$. Idem pour $\mathrm{\bf b}$ :
$$
\begin{aligned}
  \mathrm{\bf b}
  & = \sum_{I\ =\ 1}^{N_s} l(\varphi_I) \mathrm{\bf e}_I\\
  & = \sum_{p\ =\ 1}^{N_t}\sum_{I\ =\ 1}^{N_s} \underbrace{- \int_{K_p}\nabla u_\gamma^h \cdot\nabla \varphi_I}_{\text{contribution élémentaire}} \mathrm{\bf e}_I\\
\end{aligned}
$$

Afin de réduire le nombre de termes sommés, nous devons maintenant travailler localement dans chaque triangle. Pour cela, nous avons besoin d'introduire une numérotation locale de chaque sommet et utilisons une fonction $\mathrm{locToGlob}$ (*Local To Global*) permettant de basculer de l'indice local vers l'indice global telle que, pour $p \in \llbracket 1,N_t \rrbracket$ et $i \in \llbracket 1,3 \rrbracket$ :
$$\mathrm{locToGlob}\,(p,i) = I \iff \mathrm{s}_i^p = \mathrm{s}_I.$$
Nous introduisons aussi les fonctions de forme locales :
$$ \varphi_i^p = \varphi_{\mathrm{locToGlob}\,(p,i)}|_{K_p}.$$
Utilisons ces nouvelles notations en ramenant la somme sur tous les sommets du maillage à uniquement les sommets du triangle considéré :
$$
\begin{aligned}
  \mathrm{\bf A} &= \sum_{p=1}^{N_t}\sum_{i=1}^{3}\sum_{j=1}^{3} \underbrace{\int_{K_p}\nabla \varphi_j^p \cdot\nabla \varphi_i^p}_{\text{contribution élémentaire}}\ \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}^\top\mathrm{\bf e}_{\mathrm{locToGlob}\,(p,j)}\\
  \mathrm{\bf b} & = \sum_{p\ =\ 1}^{N_t}\sum_{i\ =\ 1}^{3} \underbrace{- \int_{K_p}\nabla u_\gamma^h \cdot\nabla \varphi_i^p}_{\text{contribution élémentaire}}\ \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}.
\end{aligned}
$$

Voici le pseudo-code de l'algorithme d'assemblage :

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
    b(I) += -∫_{K_p}(∇u_ɣ^h·∇ϕ_i^p)
  EndFor
EndFor
```

### Calcul des *contributions élémentaires*

Pour calculer les *contributions élémentaires*, nous devons considérer le triangle de référence $\widehat{K}$ et les fonctions de forme $\widehat{\varphi}_i \in \mathbb{P}^1(\widehat{K})$. Leur gradient est donné par :
$$
\nabla\widehat{\varphi}_1 =
\begin{pmatrix}
  -1\\
  -1
\end{pmatrix}
,
\quad
\nabla\widehat{\varphi}_2 =
\begin{pmatrix}
  1\\
  0
\end{pmatrix},
\quad
\nabla\widehat{\varphi}_3 =
\begin{pmatrix}
  0\\
  1
\end{pmatrix}.
$$
De plus, il faut considérer $\mathrm{\bf J}_p$ la matrice jacobienne de la tranformation des coordonnées de $\widehat{K}$ en celles de $K_p$ qui est donnée par :
$$
\mathrm{\bf J}_p =
\left(
\begin{array}{c c}
  x_{2}^{p} - x_{1}^{p} & x_{3}^{p} - x_{1}^{p}\\
  y_{2}^{p} - y_{1}^{p} & y_{3}^{p} - y_{1}^{p}
\end{array}
\right),
$$
avec $\mathrm{s}_i^p := (x_i^p, y_i^p)$, et dont le déterminant vaut
$$ | \det(\mathrm{\bf J}_{p})| = 2|K_p| \neq 0.$$
Les *contributions élémentaires* dans $\mathrm{\bf A}$ sont ainsi obtenus pas la relation suivante :
$$
\int_{K_p}\nabla \varphi_j^p \cdot \nabla\varphi_i^p = | K_p | (\nabla\widehat{\varphi}_j)^\top \mathrm{\bf B}_p^\top \mathrm{\bf B}_p \nabla\widehat{\varphi}_i
$$
avec $\mathrm{\bf B}_p = \left( \mathrm{\bf J}_p^\top\right)^{-1}$. Donc,
$$
\mathrm{\bf A} = \sum_{p=1}^{N_t}\sum_{i=1}^{3}\sum_{j=1}^{3} | K_p | (\nabla\widehat{\varphi}_j)^\top \mathrm{\bf B}_p^\top \mathrm{\bf B}_p \nabla\widehat{\varphi}_i\ \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}^\top\mathrm{\bf e}_{\mathrm{locToGlob}\,(p,j)}
$$
Par ailleurs, en notant,
$$
\mathcal{T}^h_{\mathcal{I}} := \{ K \in \mathcal{T}_h\ |\ K \cap (\Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}) = \empty\}
$$
l'ensemble des triangles qui ne sont pas en contact avec le bord $\Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}$, nous pouvons remarquer que, de part la définition de $u_\gamma^h$ :
$$
\mathrm{supp}\,(u_\gamma^h) \subseteq \mathcal{T}_h \backslash \mathcal{T}^h_{\mathcal{I}}.
$$
C'est-à-dire, $\forall I \in \llbracket 1, N_S \rrbracket,\forall K \in \mathcal{T}^h_{\mathcal{I}}, \mathrm{s}_I \in K \implies \mathrm{b}_I = 0$.

**TODO: heureusement gamma rad et gamma fen ne se touchent pas, sinon des cas en plus**
Si on se place sur un triangle du bord $\Gamma_{\text{Rad}} \cup \Gamma_{\text{Fen}}$, alors un triangle a soit ses trois sommets sur le bord, soit il n'en a que deux. Dans le premier cas, $u_\gamma^h$ est constante sur le triangle et donc son gradient est nul. Dans le second cas, on passe dans le triangle de référence $\widehat{K}$ où $\widehat{u_\gamma^h}$ est définie par :
$$
\widehat{u_\gamma^h} =
\left\{
\begin{aligned}
  T (\xi + \eta), & \quad \text{si $\mathrm{s}_2$ et $\mathrm{s}_3$ sont à $T$}\\
  T (1 - \xi), & \quad \text{si $\mathrm{s}_1$ et $\mathrm{s}_3$ sont à $T$}\\
  T (1 - \eta), & \quad \text{si $\mathrm{s}_1$ et $\mathrm{s}_2$ sont à $T$}
\end{aligned}
\right.
$$
où $T \in \{T_c, T_f\}$ selon le bord considéré. Les gradients sont respectivement donnés par :
$$
\nabla \widehat{u_\gamma^h} =
\begin{pmatrix}
  T\\
  T
\end{pmatrix}
,
\quad
\nabla \widehat{u_\gamma^h} =
\begin{pmatrix}
  -T\\
  0
\end{pmatrix},
\quad
\nabla \widehat{u_\gamma^h} =
\begin{pmatrix}
  0\\
  -T
\end{pmatrix}.
$$
Les *contributions élémentaires* des coefficients de $\mathrm{\bf b}$ se simplifient en :
$$
\begin{aligned}
  \int_{K_p}\nabla u_\gamma^h \cdot\nabla \varphi_i^p
  & = | \det(\mathrm{\bf J}_p) | \int_{\widehat{K}} \left( \nabla \widehat{u_\gamma^h} \right)^\top \mathrm{\bf B}_p^\top \mathrm{\bf B}_p \nabla \widehat{\varphi}_i\\
  &= | K_p | (\nabla \widehat{u_\gamma^h})^\top \mathrm{\bf B}_p^\top \mathrm{\bf B}_p \nabla\widehat{\varphi}_i

\end{aligned}
$$
Ainsi, sans tenir compte des triangles où les termes sont nuls, $\mathrm{\bf b}$ se récrit :
$$
\mathrm{\bf b} = - \sum_{p\ =\ 1}^{N_t}\sum_{i\ =\ 1}^{3} | K_p | (\nabla \widehat{u_\gamma^h})^\top \mathrm{\bf B}_p^\top \mathrm{\bf B}_p \nabla\widehat{\varphi}_i\ \mathrm{\bf e}_{\mathrm{locToGlob}\,(p,i)}.
$$

## Implantation

Utilisation des format *COO* et *CSR* pour stocker les matrices et faire les calculs.

Pas besoin de faire la matrice de masse (élémentaire) et la mettre dans l'assemblage.

Pas besoin de faire la matrice de rigidité élémentaire générique puisque formule simplifiée.

Pas besoin de faire de la quadrature.

Informatiquement, nous devons donc rendre les lignes et colonnes associées aux degrés de liberté de Dirichlet, nulles, sauf sur la diagonale avec la valeur 1. Cette opération peut être effectuée après l'assemblage de la matrice ou lors de l'algorithme directement.
Pour cela, nous parcourons les noeuds I du domaine de Dirichlet. Puis, dans
la liste des indices ligne de triplets, dès qu’une occurence à I est obtenu,
la valeur de ce triplet est mise à 0. Il ne faut pas oublier, à la fin,
d’ajouter un triplet (I,I,1) correspondant au terme diagonal et de modifier
le coefficient b[I] = g(x,y)