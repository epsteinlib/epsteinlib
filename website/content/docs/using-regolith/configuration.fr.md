---
title: "Configuration"
linkTitle: "Configuration"
weight: 4
description: >
  Effectuer des changements dans les thèmes et comportements de Regolith
next: /docs/howtos
---

## Paquets recommandés pour les nouveaux utilisateurs

### Thèmes (looks)

Installer tous les thèmes officiels de Regolith:

```bash
sudo apt install regolith-look-*
```

NOTE: Cela va télécharger plusieurs fonts et images. Si vous souhaitez minimiser l'utilisateur de votre disque, vous pouvez désinstaller les thèmes que nous n'utilisez pas.

### Indicateurs de statut

La commande suivante installe les indicateurs de statut recommandés pour les nouveaux utilisateurs:

```bash
sudo apt install i3xrocks-focused-window-name i3xrocks-rofication i3xrocks-info i3xrocks-app-launcher i3xrocks-memory
```

Pour supprimer certains de ces indicateurs plus tard, désinstallez-les simplement de votre système:

```bash
sudo apt remove i3xrocks-info
```

#### Statut des batteries pour les laptops

Pour voir le statut actuel de la batterie sur votre laptop:

```bash
sudo apt install i3xrocks-battery
```

### Compositors

Les compositors gèrent les effets visuels de vos fenêtres.
Par défaut, Regolith 2.x n'installe pas de compositor spécifique pour limiter les problèmes avec des composants graphiques (pilotes et/ou matériel).
Les utilisateurs qui souhaitent plus d'effets visuels et qui ont le matériel graphique suffisant et bien géré par Linux peuvent installer l'un des compositors suivant:

| Paquet                            | Notes                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------ |
| `regolith-compositor-picom-glx`   | Recommandé                                                                           |
| `regolith-compositor-compton-glx` | Un ancien compositor qui peut éventuellement mieux fonctionner sur certains systèmes |
| `regolith-compositor-xcompmgr`    | Un ancien compositor qui peut éventuellement mieux fonctionner sur certains systèmes |

Voir [ce guide](/docs/howtos/customize-compositor) pour plus de détails.

## Fond d'écran

La plupart des thèmes Regolith fournisse une image et une couleur par défaut pour le bureau.
Le fond d'écran peut être modifié en spécifiant le chemin vers une image via une valeur Xresources: `regolith.wallpaper.file`

Exemple, en considérant que le fichier `/usr/share/backgrounds/hardy_wallpaper_uhd.png` existe:

```bash
echo "regolith.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

Pour changer la manière dont l'image est traitée (zoom, découpage, etc..), configurer la valeur `regolith.wallpaper.options`.
Les options sont les suivantes: `none`, `wallpaper`, `centered`, `scaled`, `stretched`, `zoom`, `spanned`

Exemple:

```bash
echo "regolith.wallpaper.options: zoom" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

Pour spécifier une couleur plutôt qu'une image, utiliser la valeur `regolith.wallpaper.color.primary`:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

Pour spécifier une seconde couleur et obtenir un dégradé:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.secondary: green" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.shading.type: vertical" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Fond d'écran de verouillage

Le fond d'écran de verouillage peut être géré de la même manière que le fond d'écran du bureau.
Les valeurs Xresources sont les mêmes mais avec `lockscreen` plutôt que `wallpaper` comme préfixe (`regolith.lockscreen.wallpaper...`)
Par exemple:

```bash
echo "regolith.lockscreen.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Désactiver la gestion du fond d'écran

Si vous souhaitez gérer le fond d'écran autrement que par Regolith, fixez simplement des valeurs vides pour l'image et la couleur:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: " >> ~/.config/regolith3/Xresources
```

## Indicateurs de statut dans la barre

Les indicateurs de statut comme la charge CPU, la date et l'heure, les notifications, la météo et d'autres informations système peuvent être ajoutés ou retirés en installant des paquets.
Par exemple, pour montrer l'indicateur de statut de la batterie de votre laptop, exécutez simplement `sudo apt install i3xrocks-battery` et rafraîchir la session avec
`regolith-look refresh`.

Pour voir les indicateurs disponibles, lancez `apt search ^i3xrocks-` ou cherchez `i3xrocks-` dans votre gestionnaire de paquets graphique favoris comme [Synaptic](https://help.ubuntu.com/community/SynapticHowto).

Vous pouvez trouver [une documentation pour configurer]({{< ref "/docs/howtos/add-remove-blocklets.md" >}}) ces indicateurs de statut, aussi appelés "blocklets".

## Thèmes

Les couleurs, fond d'écran, fenêtre, barre et autres effets visuels sont regroupés dans Regolith dans des thèmes ou "looks".
Les thèmes sont une manière simple de changer entièrement l'apparence du bureau.
Comme les indicateurs de statut dans la barre, les thèmes sont des paquets qui sont installés ou retirés comme n'importe quel autre paquet.
Par convention, les paquets de thème utilise le format suivant: `regolith-look-<name>`.
`apt` ou un gestionnaire de paquet graphique peut être utilisé pour chercher les thèmes disponibles.
Pour sélectionner un thème installé, le raccourcis {{< keys "super,alt,l" >}} présente le dialogue suivant:

![](/images/v-tour/regolith-ilia-look-selector.png)

Alternativement, la commande de terminal `regolith-look` peut être utilisé pour changer de thème et rafraîchir la session active en sélectionnant un thème.
Voici un exemple qui sélectionne le thème `gruvbox`.

```console
$ apt search ^regolith-look-
[...]
regolith-look-gruvbox/unknown,now 0.4.6-1regolith amd64
[...]
$ sudo apt install regolith-look-gruvbox
regolith-look set gruvbox
regolith-look refresh
```

Pour installer l'ensemble des thèmes:

```bash
sudo apt install regolith-look-*
```

## Fonctionnalités i3

[Depuis Regolith 2.0]({{< ref "docs/reference/configurations.md#history" >}}), une grande partie de la configuration de i3 est gérable via le gestionnaire de paquet.
Regolith fourni toutes ses configurations i3 via des "configs partielles".
Ces paquets installent leurs configurations i3 dans `/usr/share/regolith/i3/config.d`.
En installant et supprimant des paquets, la configuration i3 peut être modifié selon votre préférence.
Par défaut, wuant le paquet `regolith-desktop` est installé, ces éléments de configuration sont également installés en tant que dépendances douces:

| Paquet                       | Fonction                                                       |
| ---------------------------- | -------------------------------------------------------------- |
| regolith-i3-base-launchers   | Lancer le terminal et le navigateur                            |
| regolith-i3-default-style    | Comportement des fenêtres                                      |
| regolith-i3-navigation       | Raccourcis clavier pour la navigation                          |
| regolith-i3-networkmanager   | Fonctions réseau et WiFi                                       |
| regolith-i3-next-workspace   | Utilitaires pour les espaces de travail                        |
| regolith-i3-resize           | Raccourcis clavier pour le redimensionnement                   |
| regolith-i3-session          | Raccourcis clavier pour la session                             |
| regolith-i3-unclutter        | Cache le curseur de la souris lorsque inutilisé                |
| regolith-i3-user-programs    | Lancer éventuellemetn des programmes spécifiés dans Xresources |
| regolith-i3-workspace-config | Raccourcis clavier pour les espaces de travail                 |

Les "dépendances douces" peuvent être supprimées sans causer la suppression de paquets qui dépendant d'elles.
Cela signifie que n'importe quel paquet de la liste peut être supprimé.
Les utilisateurs peuvent modifier leur configuration de manière plus stable en remplaçant une configuration partielle par leur propre version.
Cela peut se faire en copiant la configuration à remplacer dans `~/.config/regolith3/i3/config.d/` et en supprimant le fichier original dans `/usr/share/regolith/i3/config.d` via `apt`.
Ainsi, la configuration non modifiée peut toujours recevoir des mises à jour et des bug fixes sans impacter la configuration propre à l'utilisateur.

Voir [le guide pour la configuration de i3]({{< ref "docs/howtos/customize-i3-configuration.md" >}}) pour des exemples.

### Tous les paquets de configuration i3

La liste suivante contient tous les paquets de configuration pour i3 disponibles dans Regolith 2.0:

| Paquet                       | Fonction                                                               |
| ---------------------------- | ---------------------------------------------------------------------- |
| regolith-i3-base-launchers   | Lancer le terminal et le navigateur                                    |
| regolith-i3-compositor       | Intégration du Compositor (déléguée à `regolith-compositor-<variant>`) |
| regolith-i3-default-style    | Comportement des fenêtres                                              |
| regolith-i3-ftue             | Accompagnement d'un nouvel utilisateur                                 |
| regolith-i3-gaps-partial     | Configuration pour`i3-gaps` (si installé, i3-gaps remplace i3)         |
| regolith-i3-gnome            | Raccourcis clavier pour les paramètres Regolith                        |
| regolith-i3-i3xrocks         | Configuration des indicateurs de statut de la barre i3                 |
| regolith-i3-ilia             | Intégration de `ilia`                                                  |
| regolith-i3-navigation       | Raccourcis clavier pour la navigation                                  |
| regolith-i3-networkmanager   | Intégration de `network-manager`                                       |
| regolith-i3-next-workspace   | Auto-création d'espace de travail                                      |
| regolith-i3-remontoire       | Intégration de l'ancienne vue des raccourcis clavier                   |
| regolith-i3-resize           | Raccourcis clavier pour le redimensionnement                           |
| regolith-i3-rofi             | Intégration de l'ancien lanceur d'app                                  |
| regolith-i3-rofication       | Intégration de l'ancien panneau des notifications                      |
| regolith-i3-rofication-ilia  | Intégration du panneau des notifications                               |
| regolith-i3-session          | Raccourcis clavier pour la session                                     |
| regolith-i3-snapshot         | Persistance des espaces de travail i3                                  |
| regolith-i3-unclutter        | Cache le curseur de la souris lorsque inutilisé                        |
| regolith-i3-user-programs    | Lancer éventuellemetn des programmes spécifiés dans Xresources         |
| regolith-i3-workspace-config | Raccourcis clavier pour les espaces de travail                         |

## Raccourcis clavier

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

La modification la plus courante est la touche {{< keys "super" >}}.
Regolith utilise `Xresources` comme une source de vérité pour les paramètres, qui sont lu par divers composants de l'UI.
La table des valeurs `Xresources` accessibles à l'utilisateur est [accessible ailleurs]({{< ref "xresources" >}}).
Pour changer la touche par défaut {{< keys "super" >}} (touche "Win") pour la touche "Alt", ajoutez la configuration suivante au fichier `~/.config/regolith3/Xresources`:

```yaml {filename="~/.config/regolith3/Xresources"}
wm.mod: Mod1
wm.alt: Mod4
```

**Note**: GNOME utilise aussi ses propres raccourcis clavier. Quand la session Regolith est initialisée, les conflits de raccourcis avec GNOME sont supprimés
des paramètres de l'utilisateur.

![](/images/v-tour/regolith-gnome-keybindings.png)

## Gestion du système

L'app `regolith-control-center` est un outil pour configurer les langues, la date, les affichages, le réseau et d'autres paramètres.
Ouvrz le lanceur d'app avec {{< keys "super,space" >}}, tapez `paramètres` et tapez "Enter" pour lancer l'app.
Vous pouvez aussi utiliser le raccourcis {{< keys "super,c" >}} pour le lancer directement.
Cette application peut être lancée depuis un terminal avec la commande `regolith-control-center`.

![](/images/v-tour/regolith-settings-about.png)

## Plus d'infos

Pour aller plus loin, des [guides sont disponibles]({{< ref "howtos" >}}).
Vous pouvez également lire [la référence sur Xresources]({{< ref "docs/Reference/xresources.md" >}}).
Devenez un utilisateur i3 chevronné en lisant [le guide de i3](https://i3wm.org/docs/userguide.html).
