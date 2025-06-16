---
title: 'Customize i3 Configuration'
linkTitle: 'How to Customize i3 Configuration'
description: >
    How to make changes to the way i3 looks and behaves.
---

The default configurations are stored in `/usr/share/regolith/i3/config.d`. These are loaded alphabetically. Then user configurations in `~/.config/regolith3/i3/config.d` are loaded, also in alphabetical order. Regolith's default configuration is built to be customized by setting `Xresources` variables, adding user configuration, and adding or removing default configurations via `apt`. For this reason, there are several approaches to customization of i3 in Regolith that can be used separately or in combination to achieve your configuration cleanly:

-   Use `Xresources` to override variables, e.g. keybindings, strings, program names, colors, etc.
-   Add or remove Regolith default configuration files with `apt`
-   Add your own user configuration files
-   Replace Regolith's configuration wholesale

## How to set an `Xresources` variable

The following configuration example (from `/usr/share/regolith/i3/config.d/80_compositor`) can be customized to load whatever compositor you want, without editing or overriding the existing i3 configuration file, by setting `$wm.program.compositor`.

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys. The content on this page has been updated for Regolith 3.0+. If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below. For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

```text {filename="/usr/share/regolith/i3/config.d/80_compositor"}
set_from_resource $wm.program.compositor wm.program.compositor /usr/share/regolith-compositor/init
exec_always --no-startup-id $wm.program.compositor
```

To do this, simply add a line to `~/.config/regolith3/Xresources` like:

```yaml {filename="~/.config/regolith3/Xresources"}
# Use my own home-made compositor
wm.program.compositor: /usr/local/bin/my-compositor
```

Similarly, you can override other i3 options with `Xresources` based on the variable names found in `/usr/share/regolith/i3/config.d` files:

```yaml {filename="~/.config/regolith3/Xresources"}
## Gap configuration
wm.gaps.inner.size: 1

## Border configuration
wm.window.border.size: 3
wm.client.focused.color.child_border: #AAD3E9

## Workspace names
wm.workspace.01.name: 1:FOO
```

For more, please [read a more in-depth discussion on `Xresources`]({{< ref "docs/howtos/override-xres.md" >}})

## How to uninstall a default partial

If you want to wholly replace or remove settings provided by a default partial, you can use `apt` to remove the package which provides that partial.

To identify the name of the package that provides the partial you want to remove:

```bash
dpkg -S /usr/share/doc/regolith-i3-workspace-config
```

The output contains the name of the package:

```text
regolith-i3-workspace-config: /usr/share/doc/regolith-i3-workspace-config
```

Validate the full list of files provided by this package so you know what you're removing:

```bash
dpkg -L regolith-i3-workspace-config
```

Remove the partial:

```bash
sudo apt remove regolith-i3-workspace-config
```

## How to customize a default partial

{{< callout type="warning" >}}
To customize only one of Regolith's default partials, you must first
ensure that the file `~/.config/regolith3/i3/config` does **not** exist. This
overrides Regolith's defaults entirely, which is incompatible with this how-to.
{{< /callout >}}

First, ensure that a directory exists for user config partials.

```bash
mkdir -p ~/.config/regolith3/i3/config.d
```

Then, copy the default configuration into your user configuration directory:

```bash
cp /usr/share/regolith/i3/config.d/40_workspace-config ~/.config/regolith3/i3/config.d/
```

Make your desired changes:

```bash
vim ~/.config/regolith3/i3/config.d/40_workspace-config
```

Clean up the package containing the defaults:

```bash
sudo apt remove regolith-i3-workspace-config
```

Finally, **restart i3 or log out and back in**.

## Change Regolith's defaults outside of Xresources options

If you need more configuration than the above article provides please visit [Staging Config Files]({{< ref "stage-configs" >}})

## How to replace Regolith's defaults completely

If you don't want to keep any of Regolith's defaults (for example, you're an
experienced `i3` user with a complete personal config you want to use), you can
either uninstall all of Regolith's config partials, and/or create your own root
config at `~/.config/regolith3/i3/config`.
