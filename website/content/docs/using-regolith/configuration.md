---
title: 'Configuration'
linkTitle: 'Configuration'
weight: 4
description: >
    Make changes to the way Regolith looks and behaves.
next: /docs/howtos
---

## Recommended Packages for New Users

### Looks

Install all of the official Regolith looks:

```bash
sudo apt install regolith-look-*
```

Alternatively, `apt` will list available looks that may be installed:

```bash
apt list | grep regolith-look-
```

NOTE: this may pull in fonts, themes, and wallpapers. For users looking to minimize disk usage, uninstalling looks that are unused is recommended.

See below for details on how to select a look.

### Status Indicators

The following snippet installs status indicator packages recommended for new users:

```bash
sudo apt install i3xrocks-focused-window-name i3xrocks-rofication i3xrocks-info i3xrocks-app-launcher i3xrocks-memory
```

To remove any of these indicators later, simply uninstall them from your system. Example:

```bash
sudo apt remove i3xrocks-info
```

#### Laptop Battery Status

View the current battery status on a laptop:

```bash
sudo apt install i3xrocks-battery
```

### Compositors

Compositors add visual effects to the windowing environment. By default, Regolith 3 installs the `picom` compositor. Users experiencing graphical issues due to specific hardware or driver limitations may want to install one of the following compositors instead:

| Package                           | Notes                                                     |
| --------------------------------- | --------------------------------------------------------- |
| `regolith-compositor-picom-glx`   | Default                                                   |
| `regolith-compositor-compton-glx` | An older compositor that may work better on some hardware |
| `regolith-compositor-xcompmgr`    | An older compositor that may work better on some hardware |
| `regolith-compositor-none`        | A pass-through package that does not run a compositor     |

See [this How To](/docs/howtos/customize-compositor) for more details.

## Wallpaper

Most Regolith Looks provide a default desktop color or image. The wallpaper can be set by specifying the path to the wallpaper image via the Xresources key `regolith.wallpaper.file`.

Example, assuming the file `/usr/share/backgrounds/hardy_wallpaper_uhd.png` is present:

```bash
echo "regolith.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

In the above command the property mentioned is set to the `/usr/share/backgrounds/hardy_wallpaper_uhd.png` change this path to any other that needs to be used.

To change the image treatment (zoom, scale, etc.) of the wallpaper image, set `regolith.wallpaper.options`. Options are: `none`, `wallpaper`, `centered`, `scaled`, `stretched`, `zoom`, `spanned`

Example:

```bash
echo "regolith.wallpaper.options: zoom" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

To specify a color rather than an image use the Xresource key `regolith.wallpaper.color.primary`:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

To specify a secondary color and gradient:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.secondary: green" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.shading.type: vertical" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Lockscreen Wallpaper

The lockscreen wallpaper can be managed in the same way as the desktop wallpaper. The keys are the same as desktop wallpaper, but with `lockscreen` prefixed to `wallpaper`, like so: `regolith.lockscreen.wallpaper...`. For example, to specify a lockscreen image:

```bash
echo "regolith.lockscreen.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Disable Wallpaper handling

If you wish to manage wallpaper externally to Regolith, simply specify empty values for wallpaper image and color:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: " >> ~/.config/regolith3/Xresources
```

## Status Bar Indicators

Status indicators such as CPU load, date and time, notifications, weather, and other system info can be added or removed by installing packages. For example, to show a status indicator for your laptop battery, simply run `sudo apt install i3xrocks-battery` and then refresh the session using `regolith-look refresh`. To find what indicators are available, run `apt search ^i3xrocks-` or search for `i3xrocks-` in your favorite package manager GUI, such as [Synaptic](https://help.ubuntu.com/community/SynapticHowto). There is [more documentation available]({{< ref "/docs/howtos/add-remove-blocklets.md" >}}) for configuring status indicators, also called "blocklets". Here is a table of status indicators as of Regolith 2.2:

| Package                      | Status Function                  |
| ---------------------------- | -------------------------------- |
| i3xrocks-app-launcher        | Shortcut to app launch dialog    |
| i3xrocks-battery             | Battery status                   |
| i3xrocks-bluetooth           | Bluetooth status                 |
| i3xrocks-cpu-usage           | CPU utilization                  |
| i3xrocks-disk-capacity       | Disk utilization                 |
| i3xrocks-focused-window-name | Title of focused window          |
| i3xrocks-info                | Shortcut to Keybindings viewer   |
| i3xrocks-key-indicator       | Keyboard toggle status           |
| i3xrocks-keyboard-layout     | Keyboard layout                  |
| i3xrocks-media-player        | Media player status              |
| i3xrocks-memory              | Memory utilization               |
| i3xrocks-microphone          | Microphone status                |
| i3xrocks-net-traffic         | Network status                   |
| i3xrocks-next-workspace      | Shortcut to create new workspace |
| i3xrocks-nm-vpn              | VPN status                       |
| i3xrocks-openvpn             | OpenVPN status                   |
| i3xrocks-rofication          | Notification status              |
| i3xrocks-tailscale           | Tailscale status                 |
| i3xrocks-temp                | System temperature               |
| i3xrocks-time                | System time                      |
| i3xrocks-updates             | System package status            |
| i3xrocks-volume              | Audio volume status              |
| i3xrocks-weather             | Local weather                    |
| i3xrocks-wifi                | Name of connected wifi AP        |

## Looks

Colors, wallpaper, window and bar layouts and other visual factors are bundled together in Regolith and called "looks". Looks provide a simple way of changing the entire look of the desktop. Like status bar indicators, looks are packaged and are installed and removed like any other software package. By convention, Look packages use the following naming format `regolith-look-<name>`. `apt` or a GUI package manager can be used to search for available looks. To select an installed look the keybinding {{< keys "super,alt,l" >}} presents a dialog:

![](/images/v-tour/regolith-ilia-look-selector.png)

Alternatively, the terminal command `regolith-look` can be used to change looks and refresh the active session with the selected look. Here's an example that switches to the `gruvbox` look:

```console
$ apt search ^regolith-look-
[...]
regolith-look-gruvbox/unknown,now 0.4.6-1regolith amd64
[...]
$ sudo apt install regolith-look-gruvbox
regolith-look set gruvbox
regolith-look refresh
```

To simply install all available looks:

```bash
sudo apt install regolith-look-*
```

## i3 Features

[Starting with Regolith 2.0]({{< ref "docs/reference/configurations.md#history" >}}), many aspects of i3 configuration are managed via the package manager. Regolith provides all of its i3 configuration via "config partials". These packages install their i3 config partials into `/usr/share/regolith/i3/config.d`. By installing and removing packages, i3 configuration can be customized for specific preferences. By default, when the `regolith-desktop` package is installed, these configuration elements are also installed as soft dependencies:

| Package                      | Function                                                | File                                                     |
| ---------------------------- | ------------------------------------------------------- | -------------------------------------------------------- |
| regolith-i3-base-launchers   | Launch terminal and browser                             | `/usr/share/regolith/i3/config.d/15_base_launchers`      |
| regolith-i3-default-style    | Window behavior                                         | `/usr/share/regolith/i3/config.d/40_default-style`       |
| regolith-i3-navigation       | Navigation keybindings                                  | `/usr/share/regolith/i3/config.d/30_navigation`          |
| regolith-i3-networkmanager   | Network and wifi functions                              | `/usr/share/regolith/i3/config.d/88_network-manager`     |
| regolith-i3-next-workspace   | Workspace utility                                       | `/usr/share/regolith/i3/config.d/40_next-workspace`      |
| regolith-i3-resize           | Resize keybindings                                      | `/usr/share/regolith/i3/config.d/50_resize-mode`         |
| regolith-i3-session          | Session keybindings                                     | `/usr/share/regolith/i3/config.d/55_session_keybindings` |
| regolith-i3-unclutter        | Hide mouse cursor when unused                           | `/usr/share/regolith/i3/config.d/86_unclutter`           |
| regolith-i3-user-programs    | Optionally launch user programs specified in Xresources | `/usr/share/regolith/i3/config.d/90_user-programs`       |
| regolith-i3-workspace-config | Workspace keybindings                                   | `/usr/share/regolith/i3/config.d/40_workspace-config`    |

Soft dependencies can be removed without causing packages that depend upon it to be removed. This means that any of the listed packages can be removed. Users can make tweaks to configuration in a more stable manner by replacing a default configuration partial with their own version. This can be achieved by copying the partial to be customized into `~/.config/regolith3/i3/config.d/` and removing the original file in `/usr/share/regolith/i3/config.d` via `apt`. In this way, the unmodified config partials can continue to get updates and bug fixes without impacting user-specific configuration.

See the [i3 configuration customization how-to]({{< ref "docs/howtos/customize-i3-configuration.md" >}}) for examples.

### All i3 Configuration Packages

The following contains a list of all i3 configuration packages available in Regolith 2.1:

| Package                      | Function                                                              | File                                                     |
| ---------------------------- | --------------------------------------------------------------------- | -------------------------------------------------------- |
| regolith-wm-base-launchers   | Launch terminal and browser                                           | `/usr/share/regolith/i3/config.d/15_base_launchers`      |
| regolith-i3-compositor       | Compositor integration (delegates to `regolith-compositor-<variant>`) | `/usr/share/regolith/i3/config.d/80_compositor`          |
| regolith-i3-default-style    | Window behavior                                                       | `/usr/share/regolith/i3/config.d/40_default-style`       |
| regolith-wm-ftue             | First time user experience                                            | `/usr/share/regolith/i3/config.d/84_ftue`                |
| regolith-i3-gaps-partial     | `i3-gaps` configuration (if installed, i3-gaps will replace i3)       | `/usr/share/regolith/i3/config.d/35_i3-gaps`             |
| regolith-i3-gnome            | Regolith Settings keybindings                                         | `/usr/share/regolith/i3/config.d/60_config_keybindings`  |
| regolith-i3-i3xrocks         | i3 bar status indicators configuration                                | `/usr/share/regolith/i3/config.d/70_bar`                 |
| regolith-i3-ilia             | `ilia` integration                                                    | `/usr/share/regolith/i3/config.d/20_ilia`                |
| regolith-wm-navigation       | Navigation keybindings                                                | `/usr/share/regolith/i3/config.d/30_navigation`          |
| regolith-wm-networkmanager   | Network and wifi functions                                            | `/usr/share/regolith/i3/config.d/88_network-manager`     |
| regolith-i3-next-workspace   | Workspace utility                                                     | `/usr/share/regolith/i3/config.d/40_next-workspace`      |
| regolith-i3-remontoire       | Legacy keybindings view integration                                   | `/usr/share/regolith/i3/config.d/45_remontoire`          |
| regolith-wm-resize           | Resize keybindings                                                    | `/usr/share/regolith/i3/config.d/50_resize-mode`         |
| regolith-i3-rofi             | Legacy desktop launcher integration                                   | `/usr/share/regolith/i3/config.d/20_rofi`                |
| regolith-i3-rofication       | Legacy notification viewer integration                                | `/usr/share/regolith/i3/config.d/82_rofication`          |
| regolith-wm-rofication-ilia  | Notification viewer integration                                       | `/usr/share/regolith/i3/config.d/82_rofication-ilia`     |
| regolith-i3-session          | Session keybindings                                                   | `/usr/share/regolith/i3/config.d/55_session_keybindings` |
| regolith-i3-snapshot         | i3 workspace persistence                                              | `/usr/share/regolith/i3/config.d/75_i3-snapshot`         |
| regolith-wm-swap-focus       | Toggle back to last focused window                                    | `/usr/share/regolith/i3/config.d/40_i3-swap-focus`       |
| regolith-i3-unclutter        | Hides mouse when not in use                                           | `/usr/share/regolith/i3/config.d/86_unclutter`           |
| regolith-wm-user-programs    | Optionally launch user programs specified in Xresources               | `/usr/share/regolith/i3/config.d/90_user-programs`       |
| regolith-wm-workspace-config | Workspace keybindings                                                 | `/usr/share/regolith/i3/config.d/40_workspace-config`    |

## Keybindings

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys. The content on this page has been updated for Regolith 3.0+. If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below. For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

The most common keybinding change is the {{< keys "super" >}} key. Regolith uses `Xresources` as the canonical source of truth for settings, which are read by various UI components. The table of `Xresources` keys open to user configuration [is available elsewhere]({{< ref "xresources" >}}). To change the default {{< keys "super" >}} binding from the "windows" key to "alt", add the following line to the file `~/.config/regolith3/Xresources`:

```yaml {filename="~/.config/regolith3/Xresources"}
wm.mod: Mod1
wm.alt: Mod4
```

**Note**: GNOME also has its own set of keybindings. When the Regolith session is first initialized, the conflicting GNOME keybindings are removed from the user settings. GNOME keybindings can be managed in the Regolith Settings app.

![](/images/v-tour/regolith-gnome-keybindings.png)

## System Management

The `regolith-control-center` app is the tool to configure locale, date, displays, networking and various other settings. Launch it via the app launcher with {{< keys "super,space" >}}, type `settings`, and hit enter to launch the app. The direct keybinding is {{< keys "super,c" >}}. The application may also be launched from a terminal with command `regolith-control-center`.

![](/images/v-tour/regolith-settings-about.png)

## Further Reading

To dig deeper, find which [Howtos]({{< ref "howtos" >}}) are available, or read the [Xresources reference]({{< ref "docs/Reference/xresources.md" >}}). Become an i3 power user by reading their [user guide](https://i3wm.org/docs/userguide.html).

If you need more advanced configuration than Xresources provides, using i3's configuration find how to here [Staging Config Files]({{< ref "stage-configs" >}})
