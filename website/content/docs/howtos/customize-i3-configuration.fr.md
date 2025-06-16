---
title: "Personnaliser la config i3"
linkTitle: "Personnaliser la config i3"
description: >
  Comment faire des changements sur le look et les comportements de i3
---

Les configurations par défaut sont stockées dans `/usr/share/regolith/i3/config.d`.
Ces fichiers sont chargés de selon l'ordre alphabétique. Ensuite, les configurations de l'utilisateur dans `~/.config/regolith3/i3/config.d` sont chargées, également dans l'ordre alphabétique.
La configuration par défaut de Regolith est faite pour être modifiée en utilisant des variables `Xresources`, en ajoutant des configurations utilisateurs ou en ajoutant/supprimant les configurations par défaut avec `apt`.
Pour cette raison, il y a plusieurs manières de personnaliser i3 dans Regolith. Ces manières peuvent être utilisées conjointement ou séparement pour avoir une configuration propre:

- Utiliser `Xresources` pour ré-écrire les variables. Par exemple les raccourcis clavier, les labels, le nom des programmes, les couleurs,...
- Ajouter ou supprimer à la configuration par défaut de Regolith en utilisant `apt`
- Ajouter vos propres configurations utilisateur
- Remplacer la configuration Regolith de manière globale

## Comment fixer une variable `Xresources`

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

La configuration dans l'exemple suivant (récupérée de `/usr/share/regolith/i3/config.d/80_compositor`) peut être modifiée pour charger le compositor que l'on veut sans modifier ou ré-écrire la configuration existante de i3 en fixant la variable $wm.program.compositor`.

```text {filenamme="/usr/share/regolith/i3/config.d/80_compositor"}
set_from_resource $wm.program.compositor wm.program.compositor /usr/share/regolith-compositor/init
exec_always --no-startup-id $i3-wm.program.compositor
```

Pour faire cela, ajouter simplement une ligne à `~/.config/regolith3/Xresources`:

```yaml {filename="~/.config/regolith3/Xresources"}
# Utiliser mon propre programme compositor
wm.program.compositor: /usr/local/bin/my-compositor
```

De la même manière, vous pouvez ré-écrire d'autres paramètres i3 avec `Xresources` en utilisant les noms de variables trouvés dans les fichiers du dossier `/usr/share/regolith/i3/config.d`:

```yaml {filename="~/.config/regolith3/Xresources"}
## Configurer les espacements
wm.gaps.inner.size: 1

## Configurer les bordures
wm.window.border.size: 3
wm.client.focused.color.child_border: #AAD3E9

## Configurer les noms des espaces de travail
wm.workspace.01.name: 1:FOO
```

Pour plus de détails, merci de lire [la documentation concernant `Xresources`]({{< ref "docs/howtos/override-xres.md" >}}).

## Comment désinstaller une configuration partielle par défaut

Si vous souhaitez remplacer entièrement ou supprimer des paramètres fournis par une configuration spécifique, vous pouvez utiliser `apt` pour supprimer le paquet qui fourni cette configuration.

Pour identifier le nom du paquet qui donne cette configuration:

```bash
dpkg -S /usr/share/doc/regolith-i3-workspace-config
```

La sortie contient le nom du paquet:

```yaml
regolith-i3-workspace-config: /usr/share/doc/regolith-i3-workspace-config
```

Pour voir la liste complète des fichiers fournis par ce paquet (et savoir ce que l'on retire):

```bash
dpkg -L regolith-i3-workspace-config
```

Supprimer la configuration partielle:

```bash
sudo apt remove regolith-i3-workspace-config
```

## Comment personnaliser une configuration

{{< callout type="warning" >}}
Pour personnaliser une unique configuration par défaut partielle de Regolith, vous devez d'abord vous assurer
que le fichier `~/.config/regolith3/i3/config` n'existe pas.
Cela va ignorer entièrement la configuration par défaut de Regolith, ce qui n'est pas compatible avec ce guide.
{{< /callout >}}

Premièrement, assurez-vous qu'un dossier existe dans la configuration utilisateur:

```bash
mkdir -p ~/.config/regolith3/i3/config.d
```

Ensuite, copier la configuration par défaut dans le dossier de configuration de l'utilisateur:

```bash
cp /usr/share/regolith/i3/config.d/40_workspace-config ~/.config/regolith3/i3/config.d/
```

Faites les changements voulus:

```bash
vim ~/.config/regolith3/i3/config.d/40_workspace-config
```

Retirez le paquet contenant la configuration par défaut:

```bash
sudo apt remove regolith-i3-workspace-config
```

Finalement, redémarrer i3 ou sortez de la session et re-connectez vous pour appliquer les changements.

## Comment remplacer la configuration par défaut de Regolith globalement

Si vous ne souhaitez par garder la configuration de Regolith par défaut (par exemple, si vous avez des l'expérience avec `i3` avec une configuration à vous que vous souhaitez utiliser), vous pouvez désinstaller toutes les configurations de Regolith, et/ou créer votre propre configuration racine dans `~/.config/regolith3/i3/config`.
