---
title: "Modifier les Xresources"
linkTitle: "Modifier les Xresources"
description: >
  Apprenez à organiser des copies utilisateur des fichiers de configuration Regolith
---

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

Regolith repose sur le système [Xresources](https://en.wikipedia.org/wiki/X_resources) pour fournir une interface de configuration solide.
Pour change les valeurs `Xresources`, Regolith peut être modifié pour mettre à jour l'interface utilisateur, spécifier des comportements, ou définir des formats spécifiques pour certains blocs de la barre.

{{< callout type="info" >}}
Dans Regolith, les valeurs `Xresources` peuvent être changées avec deux méthodes: en fournissant vos propres fichiers `Xresources` et en ré-écrivant les valeurs `Xresources` existantes. La seconde approche est fortement recommandée car elle permet d'avoir des mises à jour de Regolith moins sensibles aux configurations changeantes.
Si vous copiez l'arbre `Xresources` entièrement, vous aurez besoin d'intégrer à la main chaque changement de configuration dans les futures versions.
{{< /callout >}}

## Initialisation

Quand une session Regolith démarre, les Xresources sont charger dans l'ordre suivant:

| Fichier                          | Optionel | Notes                                           |
| -------------------------------- | -------- | ----------------------------------------------- |
| `~/.Xresources`                  | Y        | Pour les paramètres non- Regolith               |
| `/usr/share/regolith-look/**/`   | N        | Les `Xresources` de Regolith liées à des thèmes |
| `~/.config/regolith3/Xresources` | Y        | Ré-écritures de configurations par défaut       |

Il est recommandé d'utiliser `~/.config/regolith3/Xresources` pour les modifications car cela évite d'avoir des paramètres redondants et facilite la maintenance dans le temps.

## Déterminer quelles valeurs doivent être changées

L'outil `xrdb` peut être utilisé pour lister les valeurs `Xresouces`.
Il y a [un certain nombre de valeurs existantes]({{< ref "/docs/Reference/xresources.md" >}}) dans la version R1.4.1:

```console
$ xrdb -query
[...]
gnome.icon.theme:	Moka
gnome.terminal.font:	JetBrains Mono 12
gnome.terminal.scrollbar:	never
gnome.wallpaper:	/usr/share/backgrounds/ESP_016895_1525_desktop.jpg
gnome.wm.theme:	Ayu-Mirage-Dark
[...]
```

## Exemples

Notez que les commandes ci-dessous ajoutent du texte à un fichier. Dans lancer la commande plusieurs fois va générer des lignes dupliquées.

### Exemple - Modifier l'UTI pour utiliser High DPI Screens

En utilisant le fichier `~/.config/regolith3/Xresources`, vous aurez uniquement à spécifier les valeurs que vous souhaitez changer.
L'outil `xrdb` peut être utilisé pour détemriner les valeurs actuels des variables.

1. Créer ou ajouter la valeur suivante au fichier `~/.config/regolith3/Xresources`:

```yaml {filename="~/.config/regolith3/Xresources"}
Xft.dpi: 192
```

2. Recharger la configuration `Xresources`:

```bash
regolith-look refresh
```

3. Ouvrir un nouveau terminal pour voir les changements

`192` est juste un exemple que vous pouvez adapter à vos besoins.

### Exemple - Changer la position de la barre i3

```console
$ xrdb -query | grep position
wm.bar.position:	bottom
$ echo "wm.bar.position:	top" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### Exemple - Changer le thème GTK

```console
$ xrdb -query | grep gtk
gnome.gtk.theme:	Ayu-Mirage-Dark
$ echo "gnome.gtk.theme:	Adwaita" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### Exemple - Désactiver la barre d'état

```bash
echo "wm.bar.trayoutput:	none" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Exemple - Utiliser Alt plutôt que la touche Win comme touche Super

```bash
echo "wm.mod: Mod1" >> ~/.config/regolith3/Xresources
echo "wm.alt: Mod4" >> ~/.config/regolith3/Xresources
```

Ensuite, recharger i3 pour que le changement prenne effet.

### Exemple - Lancer `nm-applet` quand i3 démarre

Certains utilisateurs préfèrent utiliser l'application `nm-applet` pour configurer et gérer le réseau sans-fil (depuis Regolith 1.5, `nm-applet` est exécuté par défaut en arrière-plan). Le fichier de configuration i3 peut être modifié pour lancer certains programmes au lancement. Mais, plutôt que de copier tout le fichier, il est possible de fixer jusqu'à trois programmes via Xresources sans avoir à changer le fichier de configuration i3. Pour cela, assurez-vous que la barre d'état système est activé (voir ci-dessus).

```bash
echo "wm.program.1: /usr/bin/nm-applet" >> ~/.config/regolith3/Xresources
```

Ce changement nécessite que vous vous déconnectiez puir reconnectiez pour prendre effet.

## Informations supplémentaires

Voir [la page de référence sur les configurations]({{< ref "/docs/Reference/configurations.md" >}}) pour plus de détails concernant la configuration de Regolith.
