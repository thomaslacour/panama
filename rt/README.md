
rt.py - coefficients de transmission et de réflexion.

Synopsis
--------

python rt.py

Dependence
----------

Les modules pythons standards suivants sont requis : - matplotlib, numpy, scipy, cmath, re, csv, pandas

Les modules non standards suivants sont requis : - mst, material

Description
-----------

`rt.py` calcule le coefficient de transmission (T) et le coefficient de réflexion (R) d'une plaque entre deux milieux, en utilisant le formalisme des matrices de transfert.

L'architecture du code de calcul est la suivante.

```console
    src/
    ├ material/
    ├ mst/
    ├ rt/
    │   ├ actions.py
    │   ├ __init__.py
    │   ├ config_file_reader.py
    │   ├ display.py
    │   ├ validation.py
    │   └ waves.py
    ├ main.py
    └ conf.txt
```

Le fichier `main.py` est le programme principal à exécuter avec python.

La configuration est spécifiée dans le fichier texte `conf.txt`.

Le code de calcul s'appuie sur l'utilisation de la base de données de matériaux "db.py" du module `material` ainsi que sur le module "mst" ("multiple scattering theory") pour le calcul des propriétés dynamiques effectives d'un milieu hétérogène.

L'exécution du code produit un fichier `data.csv` contenant les résultats de la simulation sous la forme d'un tableau.

Configuration
-------------

La configuration choisie est défini dans le fichier texte `conf.txt` se trouvant obligatoirement à la racine du répertoire. La syntaxe du fichier de configuration est la suivante :

```console
    # commentaire1
    parametre1=valeur1
    parametre2=valeur2
    # commentaire2
    parametre3=valeur3
```

Toute ligne commençant par `#` est considérée comme un commentaire et ne sera pas interprétée par le programme.

Les unités possibles pour les paramètres sont un seul choix parmi : (MHz, mm et µs) ou (kHz, m et ms) ou (Hz, km et s)

### Layers properties

Les différentes épaisseurs du multicouche sont précisées les unes à la suite des autres avec la syntaxe suivante (cas d'un milieu homogène) :

```console
   layer#1 = {nom_du_materiaux1: épaisseur 1}
   layer#2 = {nom_du_materiaux2: épaisseur 2}
   ...
   layer#N = {nom_du_materiauxN: épaisseur N}
```

Dans le cas d'un milieu hétérogène, les propriétés effectives seront calculées avec le code d'homogénéisation "mst". Le nom du matériau est obligatoirement "meta" et la syntaxe à respecter est la suivante :

```console
   layer#3 = {meta: épaisseur, Mat: matrice, Inc: inclusion, rmean: rayon, phi: fraction_volumique_en_%, poly: polydispersité_en_%}
```

Les valeurs "matrice" et "inclusion" sont des chaînes de caractères, "rayon" est un nombre.

Les matériaux disponibles par défaut sont :

> alu, eau, PU, air, acier ou Silicone_poreux.

Leurs propriétés sont définies dans le fichier `db.py`.

Les milieux extérieurs sont définis selon la syntaxe suivante (les milieux amont et aval étant respectivement les milieux gauche et droite)

```console
   halfspace_left  = nom_du_milieu_à_gauche
   halfspace_right = nom_du_milieu_à_droite
```

Les noms doivent être choisis parmi la liste des matériaux disponibles (définie dans le fichier `db.py`).

### Acoustical parameters

Deux ensembles de paramètres doivent être précisés : les paramètres fréquentiels et les paramètres pour l'angle d'incidence. Les instructions susceptibles d'être renseignées sont les suivantes :

```console
   f_min     = nombre
   f_max     = nombre
   f_fix     = nombre
   f_num     = nombre
   theta_min = nombre
   theta_max = nombre
   theta_fix = nombre
   theta_num = nombre
```

Selon le calcul demandé, tous les paramètres ne sont pas nécessaires (voir exemple). Pour calculer les valeurs des coefficients R et/ou T en incidence normale, seul "theta\_fix" est requis avec "f\_min", "f\_max" et "f\_num".

Actions
-------

Un seul choix est possible parmi différentes actions. La syntaxe est alors la suivante :

```console
   todo = action
```

La valeur "action" peut être un nombre ou une suite de plusieurs caractères. L'ensemble des actions possibles est implémenté dans le fichier "rt.py" (voir architecture). Les actions possibles sont :

- `todo=0` : calcule de R et T en fonction de la fréquence,

- `todo=1` : calcule de R et T en fonction de l'angle d'incidence,

- `todo=2` : calcule de R et T en fonction de l'angle d'incidence et de la fréquence,

- `todo=PROPN` : calcule des propriétés mécanique d'une couche particulière, où `N` désigne le numéro de la couche (layer\#N). Renvoie de la célérité de phase, de l'atténuation, de la masse volumique, de l'impédance et du module élastique en sortie.

Examples
--------

Pour calculer le coefficient de transmission et le coefficient de réflexion d'un bicouche en incidence normale en fonction de la fréquence, le fichier de configuration contiendra les instructions :

```console
   # -- begin of config file -- #
   layer#1 = {meta: 2, Mat: PU, Inc: air, rmean: 0.02, phi: 3, poly: 10}
   layer#2 = {steel: 2}
   halfspace_left  = water
   halfspace_right = water
   f_min=0
   f_max=1
   f_num=101
   theta_fix=0
   todo=0
   # -- end of config file -- #
```

Author
------

[Thomas Lacour], [Olivier Poncelet]

  [Thomas Lacour]: mailto:thomas.lacour@u-bordeaux.fr
  [Olivier Poncelet]: mailto:olivier.poncelet@u-bordeaux
