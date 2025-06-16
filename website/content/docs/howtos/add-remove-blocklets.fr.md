---
title: "Modifier les indicateurs dans la barre"
linkTitle: "Modifier les indicateurs dans la barre"
description: >
  Pour faire quelques changements dans les indicateurs de statut dans la barre
---

## Ajouter ou retirer des indicateurs de statut

### Voir les indicateurs de statut disponibles

```bash
apt search i3xrocks-
```

### Installer un nouvel indicateur de statut

```bash
sudo apt install i3xrocks-memory
regolith-look refresh
```

## Personnaliser la barre i3 {#customize-bar}

Chaque indicateur de statut dans la barre est géré par un fichier.
Le nom du fichier correspond à la position de l'indicateur dans la barre. Par exemple, un fichier avec le nom commençant par `10_` signifie que l'indicateur de statut correspondant sera placer avant un indicateur qui a un fichier de configuration commençant par `20_`.

Dans chaque fichier, on peut trouver des informations utilisées par `i3xrocks` pour executer un script et récupérer les données qui sont affichées dans la barre.
Certains des indicateurs ont des paramètres qui peuvent être ajustés à votre convenance.
Par exemple, si vous préférez que le statut de la batterie change plus vouent, vous pouvez modifier l'interval de polling.

La première chose à faire pour personnaliser la barre est de recopier les indicateurs que vous souhaitez voir à partir de `/usr/share/i3xrocks/conf.d/` vers `~/.config/regolith3/i3xrocks/conf.d`.
Chaque fichier dans `~/.config/regolith3/i3xrocks/conf.d`. peut être modifier selon vos besoins.
Pour changer l'ordre des indicateurs, changez simplement le nom des fichiers.

Une fois que vous avez fait des changements, il est nécessaire de rafraîchir la session pour que la barre se mette à jour en fonction de la nouvelle configuration.
Par exemple, pour changer l'ordre des blocs de la batterie et du trafic réseau et cacher les notifications, suivez ces étapes:

```console
$ ls /usr/share/i3xrocks/conf.d/
01_setup
30_net-traffic
80_battery
80_rofication
90_time
$ mkdir -p ~/.config/regolith3/i3xrocks/conf.d
$ cd /usr/share/i3xrocks/conf.d/
$ cp 01_setup ~/.config/regolith3/i3xrocks/conf.d/01_setup
$ cp 80_battery ~/.config/regolith3/i3xrocks/conf.d/30_battery
$ cp 30_net-traffic ~/.config/regolith3/i3xrocks/conf.d/80_net-traffic
$ regolith-look refresh
```

**Note**: Si un bloc existe dans le dossier utilisateur `~/.config/regolith3/i3xrocks/conf.d`, la configuration par défaut dans `/usr/share/i3xrocks/conf.d/` sera ignorée.

## Informations complémentaires

La [page de référence des configuration]({{< ref "/docs/Reference/configurations.md" >}}) donne plus de détails à propres des fichiers de configuration utilisés par Regolith Linux.
