---
title: "Configurations"
linkTitle: "Configurations"
weight: 2
description: >
  Information générale concernant les fichiers de configuration utilisés par Regolith
---

## Référence des fichiers de configuration

| **Composant**              | **Configuration par défaut**                                 | **Configuration utilisateur**         | **Notes**                                                                                                                                                                                                                         |
| -------------------------- | ------------------------------------------------------------ | ------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| i3 / i3-gaps               | `/etc/regolith/i3/config`, `/usr/share/regolith/i3/config.d` | `~/.config/regolith3/i3`              | _Étendre_ la configuration Regolith2 en créant `config.d` dans le dossier de configuration de l'utilisateur ou globalement en créant un dossier `config`. [Plus d'infos]({{< ref "docs/howtos/customize-i3-configuration.md" >}}) |
| Xresources                 | `/usr/share/regolith-look/default`                           | `~/.config/regolith3/Xresources`      | `~/.Xresources` est aussi chargé mais utilisé pour des configurations qui peuvent être nécessaire à d'autres sessions de bureau. [Plus d'infos]({{< ref "docs/howtos/override-xres.md" >}})                                       |
| Bar status indicators      | `/usr/share/i3xrocks/conf.d/`                                | `~/.config/regolith3/i3xrocks/conf.d` | Chaque bloc a son propre fichier. Les configurations sont chargées dans l'ordre alphabétique de leur nom de fichier. [Plus d'infos]({{< ref "docs/howtos/add-remove-blocklets.md" >}})                                            |
| Picom (Compton) Compositor | `/etc/regolith/picom/config`                                 | `~/.config/regolith3/picom/config`    | [Plus d'infos]({{< ref "docs/howtos/customize-compositor.md" >}})                                                                                                                                                                 |

{{< callout type="warning" >}}
Ajouter une config i3 dans `~/.config/regolith3/i3/config` va écraser toutes les configurations Regolith, y compris pour les Xresources et les thèmes.
{{< /callout >}}

## Historique

### Configurations partielles

Regolith 2.0 ajoute la gestion de configurations partielles suite à l'ajout de la fonctionnalité dans i3 récemment ([i3 4.20](https://i3wm.org/downloads/RELEASE-NOTES-4.20.txt)). Cela est réalisé en utilisant les [inclusions](https://i3wm.org/docs/userguide.html#include).
La gestion des inclusions permet à Regolith 2.0 d'utiliser le gestionnaire de paquet pour ajouter ou supprimer des configurations partielles par défaut.
En installant ou en supprimant des paquers, la configuration de i3 peut être personnalisée selon les préférences tout en suivant les mises à jour de configuration.

Le fichier de configuration Regolith racine (`/etc/regolith/i3/config`) charge les configurations partielles de la manière suivante:

```text {filename="/etc/regolith/i3/config"}
# Inclus toutes les configuration de Regolith pour i3
include /usr/share/regolith/i3/config.d/*

# Inclus toutes les configurations de l'utilisateur pour i3
include $HOME/.config/regolith2/i3/config.d/*
```
