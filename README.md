## TODO

* [X] externaliser le calcul d'homogénéisation (calcul sur toutes les fréquences)
* [X] gerer les entrées de calculs "stacks"[frequences][couches,propriétés]
      - temps de calcul gagné (facteur 7 à 10)
* [X] généraliser sur plusieurs couches
* [X] vérifier la simulation bicouche, les BFs ne correspondent pas à une situation eau|mat|air
      - !erreur trouvée sur la taille des inclusions
- [ ] correction bug formule analytique angle dans `validation.py`. (la simu est bonne mais pas la formule analytique)
* [X] intégrer fonction "géométrie" qui génère les configurations du multicouche
- [ ] supprimer les fonctions inutiles de `rt/display.py` et revoir les fonctions d'affichage
- [X] revoir la gestion du vecteur fréquence dans le module MST (problème de shape)
- [ ] intégrer multicouche dans main.py
- [ ] externaliser et dissocier certains calcul dans `waves.py`


