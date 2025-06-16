---
title: "Config Files"
linkTitle: "Config Files"
weight: 2
description: >
  General information about the config files that Regolith uses
---

## Config file quick reference

| **Component**              | **Default Config**                                           | **User Config**                       | **Notes**                                                                                                                                                                                             |
| -------------------------- | ------------------------------------------------------------ | ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| i3 / i3-gaps               | `/etc/regolith/i3/config`, `/usr/share/regolith/i3/config.d` | `~/.config/regolith3/i3`              | _Extend_ Regolith2 config by creating `config.d` and partials in the user config dir, or _globally override_ it by creating `config`. [More]({{< ref "docs/howtos/customize-i3-configuration.md" >}}) |
| Xresources                 | `/usr/share/regolith-look/default`                           | `~/.config/regolith3/Xresources`      | `~/.Xresources` is also loaded but intended for properties that may also be required in other desktop sessions. [More]({{< ref "docs/howtos/override-xres.md" >}})                                    |
| Bar status indicators      | `/usr/share/i3xrocks/conf.d/`                                | `~/.config/regolith3/i3xrocks/conf.d` | Each block has it's own file. Configs loaded in filename alphabetical order. [More]({{< ref "docs/howtos/add-remove-blocklets.md" >}})                                                                |
| Picom (Compton) Compositor | `/etc/regolith/picom/config`                                 | `~/.config/regolith3/picom/config`    | [More]({{< ref "docs/howtos/customize-compositor.md" >}})                                                                                                                                             |

{{< callout type="warning" >}}
Adding an i3 config file in `~/.config/regolith3/i3/config` will essentially override all Regolith configuration including Xresources and looks.
{{< /callout >}}

## History

### Config partials

Regolith 2.0 adds support for config partials, as this feature came to `i3` recently in [4.20](https://i3wm.org/downloads/RELEASE-NOTES-4.20.txt). This is done by [includes](https://i3wm.org/docs/userguide.html#include). This support for includes enables Regolith 2.0's use of the package manager to add or remove default configuration partials. By installing and removing packages, i3 configuration can be customized for specific preferences while still allowing to track upstream changes for aspects of the configuration that need not vary.

The root Regolith configuration file at `/etc/regolith/i3/config` loads partials like so:

```text {filename="/etc/regolith/i3/config"}
# Include any regolith i3 partials
include /usr/share/regolith/i3/config.d/*

# Include any user i3 partials
include $HOME/.config/regolith2/i3/config.d/*
```
