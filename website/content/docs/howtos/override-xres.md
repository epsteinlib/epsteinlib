---
title: "Override Xresources"
linkTitle: "Override Xresources"
description: >
  Learn how to stage user copies of Regolith configuration files
---

Regolith relies on the [Xresources system](https://en.wikipedia.org/wiki/X_resources) to provide a consolidated interface configuration. By changing `Xresources` values, Regolith can be customized in ways such as updating the user interface, specifying custom behaviors, or defining a specific format for some of the blocklets on the bar.

{{< callout type="info" >}}
In Regolith, `Xresources` values can be changed via two methods: supplying your own `Xresources` files, and overriding the existing `Xresources` values. The latter approach is strongly recommended, as it makes upgrading to newer versions of Regolith less likely to break your customizations. If you copy the `Xresources` tree completely you'll need to integrate any breaking changes in future versions manually.
{{< /callout >}}

## Initialization

When a Regolith session is started, Xresources are loaded in the following order:

| File                            | Optional | Notes                                                                             |
| ------------------------------- | -------- | --------------------------------------------------------------------------------- |
| `~/.Xresources`                 | Y        | Intended for non-Regolith settings                                                |
| `/usr/share/regolith-look/**/`     | N        | The Regolith `Xresources` associated with looks |
| `~/.config/regolith3/Xresources` | Y        | Applies specific overrides to `Xresources` defaults                               |

It is recommended to use `~/.config/regolith3/Xresources` for customization as it doesn't require the specification of redundant settings and is easier to maintain over time.

## Determining which values can be changed

The `xrdb` tool can be used to list the existing `Xresources` values. [There is a number of existing overrides]({{< ref "/docs/Reference/xresources.md" >}}) in the R2.1 release:

```console
$ xrdb -query
[...]
gnome.terminal.scrollbar:	never
gnome.terminal.use-transparent-background:	true
gtk.document_font_name:	Sans 12
gtk.font_name:	Bitstream Vera Sans 12
gtk.icon_theme_name:	Papirus-Dark
gtk.monospace_font_name:	BitstreamVeraSansMono Nerd Font 13
gtk.theme_name:	Ayu-Mirage-Dark
wm.bar.background.color:	#1F2430
[...]
```

## Examples

Note that the commands presented below append text to a file, so running the command more than once will result in duplicate entries.

### Example - Update the UI for High DPI Screens

By using the `~/.config/regolith3/Xresources` override file, we will only need to specify the values we wish to change. The `xrdb` tool can be used to determine what current values are set to.

1. Create or add the following value to your `~/.config/regolith3/Xresources` file:

```yaml {filename="~/.config/regolith3/Xresources"}
Xft.dpi: 192
```

2. Reload the `Xresources` configuration:

```bash
regolith-look refresh
```

3. Open a new terminal to see the change take effect.

`192` is just an example value, please adjust as needed.

### Example - Change i3 bar position

```console
$ xrdb -query | grep position
wm.bar.position:	bottom
$ echo "wm.bar.position:	top" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### Example - Change GTK Theme

```console
$ xrdb -query | grep gtk
gtk.theme_name:	Ayu-Mirage-Dark
$ echo "gnome.gtk.theme:	Adwaita" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### Example - Disable the System Tray

```bash
echo "wm.bar.trayoutput:	none" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Example - Use Alt instead of Win as Super

```bash
echo "wm.mod: Mod1" >> ~/.config/regolith3/Xresources
echo "wm.alt: Mod4" >> ~/.config/regolith3/Xresources
```

Then Reload i3 for the change to take effect.  Some settings may require logging back into the session, for example anything that i3 launches as a separate process.

## Further Reading

See the [reference page for configurations]({{< ref "/docs/Reference/configurations.md" >}}) for more details about config files in Regolith.
