---
title: "1.3.1 Release Notes"
linkTitle: "1.3.1 Release Notes"
weight: 98
description: >
  Release notes for Regolith 1.3.1.
---

R1.3.1 is a quick followup to R1.3 that adds polish and usability features. Headlining features include a more concise way of overriding Xresource defaults and several enhancements to the Remontoire keybinding viewer.

{{< callout type="info" >}}
If you maintain your own i3 config file, note that the copy shipped in Regolith 1.3.1 has changed in that the comment format has been modified to be readable by Remontoire ([example](https://github.com/regolith-linux/regolith-i3-gaps-config/blob/master/config#L37)). Without updating your copy i3 config comments, you will not see any keybindings.
{{< /callout >}}

## Features

| Feature                                                  | Screenshot                                                                               | Description                                                                                                                                                                                                                                                                                                                                                   |
| -------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Remontoire 1.1.1                                         | {{< img "images/releases/regolith-remontoire-screenshot.png" "New Remontoire release">}} | Provides a live view of active i3 keybindings, introduces a machine-readable comment format for i3 config files, and fixes several layout and presentation bugs reported by users.                                                                                                                                                                            |
| Ayu Looks                                                | {{< img "images/releases/regolith-screenshot-ayu-mirage.png" "New theme" >}}             | The community contributed Ayu themes have been fully integrated into Regolith with updates to package dependencies and configurations.                                                                                                                                                                                                                        |
| Rofi polishing                                           | {{< img "images/releases/regolith-screenshot-rofi-window.png" "Rofi polishing" >}}       | Numerous updates to the Rofi themes for all provided looks.                                                                                                                                                                                                                                                                                                   |
| Easier Customization                                     |                                                                                          | Now, overriding specific `Xresources` values can be easily done by supplying only the overrides in `~/.config/regolith3/Xresources` without needing to copying and update any existing files. Additionally, any `Xresources` files in `~/.config/regolith3/Xresources.d` can be referenced without absolute path, enabling easier multi-machine configurations. |
| Utilize GNOME and `update-alternatives` for Default Apps |                                                                                          | GNOME settings for default browser and `update-alternatives` for preferred terminal are now used by Regolith. No custom configuration changes are required.                                                                                                                                                                                                   |
| `Xresources` updates                                     |                                                                                          | More UI configuration pulled out of component configs and into Xresources provides more flexibility in how custom looks can be.                                                                                                                                                                                                                               |
| Locale used for date/time format                         |                                                                                          | For those not using the US date and time formats, the bar clock will now automatically use the GNOME hour format and the locale `LC_TIME` setting to automatically display date and time in the regional form without configuration.                                                                                                                          |
| UI tweaks and polishing                                  | {{< img "images/releases/nord-dark.png" "Nord Dark theme" >}}                            | Many small updates across all UI components to further polish the interface, including a new optional dark background adapted from work by [/u/InfernusZ9](https://www.reddit.com/r/unixporn/comments/ehdj6b/i3gaps_nordic_beauty/).                                                                                                                          |

## Fixes

- Debian package dependencies for various Regolith Looks have been fixed.
- Layout issues with Remontoire have been addressed.
- Compton configuration cleaned up
