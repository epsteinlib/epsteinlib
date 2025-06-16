---
title: "Les bases"
linkTitle: "Les bases"
weight: 3
description: >
  Commencez par apprendre les bases importantes de Regolith
---

Bien que nous couvrirons l'ensemble du fonctionnement de Regolith dans d'autres parties de la documentation, cette page couvre les cas d'utilisations basiques.

## Lancer une application

Après la première connexion, l'écran contiendra une barre sur le bas et quelques raccourcis clavier dans un dialogue.
Le but de ce dialogue est d'aider l'utilisateur à trouver et se souvenir les raccourcis clavier pour diverses actions sur l'environnement de bureau.
Il peut être caché en pressant la touche `esc` et être ré-affiché avec {{< keys "super,shift,?" >}}.
Pour lancer une application, pressez {{< keys "super,shift,?" >}}. et un dialogue similaire affichera les applications installées sur le système.
Si vous commencez à taper le nom de l'application que vous souhaitez lancer, la liste sera filtrée. Vous pouvez aussi scroller dans la liste avec les touches de votre clavier.

![](/images/v-tour/regolith-ilia-apps.png)

## Naviguer vers une fenêtre

Pour se focaliser sur une application en cours d'exécution, le raccourcis {{< keys "super,ctrl,space" >}} peut être utilisé pour afficher un dialogue avec toutes les fenêtres ouvertes. Sélectionnez simplement la fenêtre voulue et vous serez déplacés vers l'espace de travail contenant cette fenêtre.

![](/images/v-tour/regolith-ilia-windows.png)

## Créer un nouvel espace de travail

Les espaces de travail (aussi appelés "bureaux virtuels" sur certains environnements de bureau) sont la première manière de gérer les fenêtres.
Chaque espace de travail correspond à un écran avec un ensemble indépendant de fenêtre. Pour créer ou naviguer vers un espace de travail, spécifiez son index avec les touches {{< keys "super,0-9" >}}.

## Trouver et lancer un fichier

Rechercher un fichier spécifique sur votre système peut être fait de différentes manières.
Regolith fourni un dialogue de recherche de fichier avec {{< keys "super,alt,space" >}} pour rapidement chercher et éventuellement afficher un fichier avec son application par défaut. Pressez 'enter' sur un élément sélectionné lancera une tentative d'ouverte avec la commande `xdg-open`.

![](/images/v-tour/regolith-ilia-files.png)

## Terminal

Le terminal est une application spéciale qui peut être lancé avec le raccourcis {{< keys "super,enter" >}}.

## Navigateur Web

Le navigateur préféré est aussi spécial et peut être lancé avec {{< keys "super,shift,enter" >}}.

## Quitter une application

{{< keys "super,shift,q" >}} fermera l'application focalisée.
Pour les applications qui ne réagissent plus, {{< keys "super,alt,q" >}} les forcera à quitter (kill) mais notez que des données non sauvées peuvent être perdues.

## Naviguer entre les espaces de travail

Le moyen principal de gérer comment les applications occupent l'écran est d'utiliser les espaces de travail.
Un espace de travail est un groupe de 0 ou plus applications qui sont toutes présentes sur le même affichage virtuel.
Les utilisateurs peuvent naviguer entre les espaces de travail pour gagner plus de place.
Les fenêtes peuvent être déplacées entre les espaces de travail.
Les espaces de travail peuvent être attachés à un écran physique ou être déplacés vers un autre écran au besoin.
Pour naviguer entre les espaces, tenez la touche {{< keys "super" >}} et appuyez sur la touche d'un numéro.
Cela va déplacer la vue vers cet espace de travail. S'il y a des fenêtres sur cette espace, elles seront listées dans la barre en bas.
Ainsi, en regardant la barre du bas, l'utilisateur peut savoir quelles sont les fenêtes actives sur l'espace de travail focalisé.

Il est également possible de se déplacer "d'avant en arrière" entre les espaces de travail en utilisant {{< keys "super,tab" >}} and {{< keys "super,shift,tab" >}}.

## Configuration du système

L'application Paramètres (aussi connu sous le nom `regolith-control-center`) est l'interface principale pour les tâches de configuration courantes comme la configuration du WiFi et du Bluetooth, la configuration des apps associée et la configuration des utilisateurs.
Cette app peut être lancée avec le raccourcis {{< keys "super,c" >}} ou depuis le lanceur d'applications en utilisant le nom `Regolith Settings`.

![](/images/v-tour/regolith-gnome-settings.png)

## Liste des raccourcis clavier

Tous les raccourcis clavier utilisables sur Regolith Desktop peuvent être vus, triés et exécutiées par `ilia`.
Ce dialogue est lançable avec {{< keys "super,shift,?" >}}.

![](/images/v-tour/regolith-ilia-keybindings.png)

## Gestion des notifications

Le système de notification de Regolith est fait pour prévenir les pop-ups qui peuvent distraire l'utilisateur de son travail.
Les notification peuvent être vues et gérées avec {{< keys "super,n" >}}.
Les raccourcis clavier pour supprimer une ou plusieurs notifications sont trouvables dans l'ongler `aide` du dialogue.
Pour les utilisateurs qui souhaitent plus de visibilité sur les notifications en cours, l'indicateur de status `i3xrocks-rofication` peut être installé via le système de gestion de paquet. Cela affichera le nombre de notifications critiques et non-vues dans la barre.
