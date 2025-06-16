---
title: "Modifier le Compositor"
linkTitle: "Modifier le Compositor"
description: >
  Comment modifier le compositor utilisé par Regolith
---

{{< callout type="info" >}}
Portez une attention particulière à cette section si vous rencontrez des problèmes visuels étranges,
des problèmes ou des performances graphiques lentes.
{{< /callout >}}

Un _compositor_ est un composant d'interface utilisateur qui applique des effets visuels
aux fenêtres avant qu'elles ne soient affichées sur l'écran.
Beaucoup d'environnements de bureau intègrent un compositor directement dans leur gestionnaire de fenêtre,
rendant cela difficile à désactiver ou à modifier.

Dans Regolith, le compositor est défini comme une extension "pluggable" avec le système de paquets.
Cela signifie qu'il est possible de changer de compositor en installant un paquet qui en contient un.
Le système de gestion de paquets s'assure qu'il n'y a pas de conflit et que toutes les dépendances sont satisfaites.

{{< callout type="warning" >}}
Dans Regolith 2.x, le compositor par défaut est "pas de compositor".
Dans Regolith 1.x, le compositor par défaut est Picom.
{{< /callout >}}

## Lister les compositors disponibles

La commande suivante liste tous les compositors configurés qui fonctionnent avec Regolith:

```bash
apt search regolith-compositor-
```

Vous devriez trouver, au minimum:

- `regolith-compositor-none` **[par défaut]**: Aucun compositor. Les meilleures performances, pas d'effet visuel.
- `regolith-compositor-picom-glx`: Utilise `picom`. Recommandé pour la plupart des utilisateurs pour résoudre les problèmes de "déchirure" d'écran et ajouter des effets visuels. Moderne et bien maintenu, fork de `compton`. Il était le compositor par défaut dans Regolith 1.6.
- `regolith-compositor-compton-glx`: Un ancêtre de `picom` qui n'est plus maintenu aujourd'hui.
- `regolith-compositor-xcompmgr`: Un compositor récent qui peut offrir des performances meilleures si les compositors complexes posent soucis.

## Installer un compositor

{{< callout type="info" >}}
À cause de la manière dont le compositor est géré par i3-wm, vous devez vous déconnectez
et vous reconnectez pour que la modification de compositor prenne effet.
Redémarrer i3 n'est pas suffisant.
sufficient.
{{< /callout >}}

Lancez la commande d'installation suivante pour changer le compositor.
Cela installera `picom`:

```bash
apt install regolith-compositor-picom
```

**Déconnectez puis reconnectez vous** et le nouveau compositor devrait être chargé.
Vous pouvez vérifier cela dans la liste des processus en cherchant `picon`:

```bash
ps aux | grep picom
```

## Modifier la configuration par défaut de Regolith

### Compton/Picom

Pour fournir votre propre configuration pour Compton/Picom, copiez la config par défaut
et créez la votre dans `~/.config/regolith3/picom/config`.
À la prochain session, la nouvelle config sera chargé à la place de celle par défaut (`/etc/regolith/picom/config`).
