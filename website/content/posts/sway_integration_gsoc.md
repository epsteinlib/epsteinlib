---
authors:
  - name: Soumya Ranjan Patnaik
    link: https://github.com/SoumyaRanjanPatnaik
    image: https://avatars.githubusercontent.com/u/75623505
date: 2022-09-25
linktitle: "Foundations for Regolith on Wayland (GSoC'22)"
title: "Foundations for Regolith on Wayland (GSoC'22)"
---

## Introduction

Regolith is a modern desktop environment designed to let you work faster by reducing unnecessary clutter and ceremony. It offers a keyboard based workflow and a tiling window manager that allows efficient management of windows, monitors and workspaces. Built on top of Ubuntu, GNOME Flashback, and i3, Regolith stands on a well-supported and consistent foundation.

In its current implementation, the use i3 window manager binds Regolith to X11. While X11 is fairly mature, it was never built for modern systems and is based on assumptions that are valid only for computers of the ’90s. Some of the issues present in X11 based window managers are low framerates, poor support for HiDPi panels, higher latency, screen tearing etc. Wayland is a display protocol that is written from scratch and aims at simplicity and extensibility Wayland has been in development for quite some time now and its implementations are now stable enough to be considered for daily usage. In fact, many desktop environments like gnome are now providing their wayland implementations as default. Sway window manager was introduced as a drop-in replacement for i3. It is based on wayland and is devoid of most of the issues present in X11 based window managers and compositors. 

>**Note:** The work described in this post was done as part of my [Google Summer of Code (GSoC) 2022](https://summerofcode.withgoogle.com/) project under the organisation [CCExtractor](https://ccextractor.org/).

## Objective

The project intends to achieve integration of swaywm with the regolith ecosystem. This means, compatibility with GNOME services, integration of a configuration management system, theming using regolith look and support for other extras provided by regolith in a wayland-based environment.

## Integrating sway with regolith

In this section, I will described everything that I did to achieve sway integration with regolith. This required modifying several existing packages as well as writing code for new scripts and binaries.

### Porting the regolith app launcher (ilia) to wayland

After having a discussion with Ken, we decided that the most important part of a tiling window manager setup from a user's perspective is the app-launcher. Regolith on X11 ships with its own app launcher called `ilia`. It also provides several features apart from launching apps. This includes showing keybindings, launching terminal commands, list applications and workspaces and search files just to name a few.

Ilia is written in `vala` and `GTK`. This means, ilia had support for wayland right out of the box. Despite this inherent compatibility, ilia seemed to crash everytime. After looking into the logs, I realised that the issue was not with `ilia` but with the way `wayland roots` (wlr) handles windows that need to be ignored by the compositor. Turns out, wlr provides an extension to the wayland protocols called `wlr-layer-shell` in order to handle panels and other applications that need to be positioned automatically. A client implementation of this protocol for GTK based applications is provided by the library `gtk-layer-shell`. To maintain support for both `X11` and `wayland`, I check the `XDG_SESSION_TYPE` environment variable and use gtk-layer-shell library only if the session type is wayland.

Currently, the wayland compatible version of `ilia` is not packaged and exists in a [separate branch](https://github.com/regolith-linux/ilia/tree/ubuntu/jammy-wayland).

### Trawl - A backend agnostic configuration manager

`Wayland` lacks support for `Xresources` natively and requires the use `XWayland`. This might not be desirable as it forces users to enable `Xwayland` forcing them to be a part of a legacy ecosystem. This is the reason why we decided to create a display backend agnostic configuration manager that can be used by various applications. To maintain compatibility with the existing regolith ecosystem it was crutial to replicate the format and functionaility XResources as much as possible. `Trawl` is a file based configuration management system written in rust that tries to fulfill these requirements.

[`Trawl`](https://github.com/regolith-linux/trawl) is comprised of two parts - Configuration Daemon (`trawld`) and Trawl Client (`trawldb`). For IPC, it uses the `DBus` specification, which allows other applications to access `trawl` resources with ease. In fact, trawl comes with libraries for [C](https://github.com/regolith-linux/libtrawldb) and [rust](https://crates.io/crates/trawldb) that make interaction with trawl clients hassle free. Trawl also tries to replicate most of the features provided by `Xresources` by replicating the functionality available using the `xrdb` CLI.

### Sway - adding the set_from_resource command

`Sway` aims to be a drop in alternative for `i3`. But there are some features which make sense only within the `X11` ecosystem, and one such feature is the `set_from_resource` command in the configs. This command allows users to set variables within the i3 configurations by querying `XResources` for specific keys. Due to the lack of `XResources` alternatives that are display backend agnostic, the creators of sway decided [not to implement it](https://github.com/swaywm/sway/issues/2336#issuecomment-469041201), even if its at the cost of breaking compatibility with i3 configurations.

Regolith relies heavily on `XResources` based configurations. This helps acheive functional configurations to some extent, such that you can replicate your window manager setup just by copying over your `XResources` and `config partials`. In order to maintain a consistent experience across the board, it was absolutely crutial to re-implement `set_from_resource` for sway.

I created a fork of sway called [`sway-regolith`](https://github.com/regolith-linux/sway-regolith) which adds the [`set_from_resource`](https://github.com/regolith-linux/sway-regolith/blob/7d4bca923950db19bc3341d960ed5862339eb1b7/sway/commands.c#L87) command. In order to communicate with the trawl demon, I used the wrapper functions provided by the [`libtrawldb`](https://github.com/regolith-linux/libtrawldb) library, which provides several wrapper functions for the DBus methods. Finally, I renamed the `sway` binary to `sway-regolith` in order to allow the fork and the upstream versions of sway to co-exist.

### Regolith-displayd

In regolith, the de facto way for managing and changing settings is the use of `gnome-control-center`. Since the displays are managed by `X11` in both `i3` and `Gnome on X11`, the display panel of `gnome-control-center` works right out of the box. However, in the case of wayland based compositors and window managers like `sway` and `mutter` on wayland, the diplays are managed by the compositor itself. So `gnome-control-center` can manage displays only when mutter is the compositor in use.

On looking into `gnome-control-center`'s codebase, I discovered that the display settings are applied by communicating with the `org.mutter.DisplayConfig` DBus interface. So I created [`regolith-displayd`](https://github.com/regolith-linux/regolith-displayd), which creates a `DBus` interface with the same name. It primarily provides two DBus methods - `GetCurrentState` and `ApplyMonitorsConfig`. The `GetCurrentState` method gets the current configuration of monitors and maps it into a format understandable by `gnome-control-center`. Similarly, calls to the `ApplyMonitorsConfig` method applies the display settings changes made from `gnome-control-center` by mapping them into a format understandable by the IPC interface of sway.

`regolith-displayd` also polls the sway IPC interface to check for changes in display configuration. If a change is detected, the `MonitorsChanged` signal is emmitted. `gnome-control-center` reacts to this message by making a call to `GetCurrentState`.

In order to achieve display configuration persistance, `regolith-displayd` generates `kanshi` profiles based on the displays that are currently in use. These configurations are then automatically applied whenever the same monitors are connected.

### Regolith-session

`regolith-session` package is probably the most important part of regolith as it essentially ties i3 with GNOME flashback. This package is also partly responsible for the `regolith-look` related functionality. The code for wayland version of the package can be found in the [`wayland-dev`](https://github.com/regolith-linux/regolith-session/tree/wayland-dev) branch.

Migrating regolith-session to wayland required work on the following sections.

#### Session initialization script (`regolith-session-wayland`)

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/bin/regolith-session-wayland) script is responsible for loading the `trawl` resources set by the user, loading current look information, reformatting pango markup for workspace names, registering the session with `gnome session manager` and launching sway with the regolith configuration file.

#### `/usr/share/applications/regolith-wayland.desktop`

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/share/applications/regolith-wayland.desktop) application file executes the `regolith-session-wayland` script. This application file will be read by `gnome-session-binary` while loading the Required Components from the session file.

#### `/usr/share/gnome-session/sessions/regolith-wayland.session`

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/share/gnome-session/sessions/regolith-wayland.session) file lists out the components required by the session in order to work correctly. If all the required components cannot be loaded, the session ends abruptly with an error. Most of these components are `gnome-settings-daemon` plugins. If the session is loaded using `gnome-session --builtin`, it launches the required components using their respective application files instead of launching their systemd services. The `regolith-sway` desktop file is also included in the required components section which launches ultimately starts the regolith session.

#### `/usr/share/wayland-sessions/regolith-sway.dektop`

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/share/wayland-sessions/regolith-sway.desktop) file is read by the display manager / login manager and shows the entry `Regolith on Wayland` in the login screen. The command `/usr/bin/gnome-session --builtin --session=regolith-wayland` is executed as soon as you login. This command looks for `/usr/share/gnome-session/sessions/regolith-wayland.session` and loads the required components.

#### `/usr/lib/regolith/regolith-session-common.sh`

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/lib/regolith/regolith-session-common.sh) provides functionality that is used in multiple scripts (related to both `X11` and `wayland` based sessions). It also contains some variables that are required by many of the scripts related to `regolith-session`. I've added several functions and environment variables that are required to provide functionality equivalent to what is provided by Regolith on X11 sessions.

#### Regolith look script (`/usr/bin/regolith-look`)

[This](https://github.com/regolith-linux/regolith-session/blob/c92813f5c836443ab0d85eb4e925a3aa2bbbfd2e/usr/bin/regolith-look) provides theming related functionallity, which is one of the defining features of the regolith ecosystem. I've added a check to determine if the current session type (X11 or Wayland). If the current session is wayland, `trawlcat` is used instead of `xrescat` and `swaymsg reload` is used instead of `i3-msg refresh`. When loading a look, if the current session is wayland, the configurations are loaded into `trawl` and if the session is X11, the resources are loaded using `xrdb` as XResources.

### Regolith Look

The functionality of the `regolith-look` script (located at `/usr/bin/regolith-look`) has been discussed in the [previous section](#regolith-look-script-usrbinregolith-look). In order to apply the theme, a loader script is required. This can either be provided by the look or the user. If no loader is provided, the [default loader script](https://github.com/regolith-linux/regolith-look-default/blob/1d50b0797c579b239664a8a96aaf7de399f0f2da/usr/share/regolith-look/default_loader.sh) (present in the `regolith-look-default` package) is used. To make this loader script compatible with `trawl` in wayland, I had to check if the current session is wayland and set the `RESOURCE_GETTER` variable to `trawlcat` accordingly.

Regolith also uses gnome-terminal by default. So in order to provide terminal theming in wayland, the [`regolith-init-term-profile`](https://github.com/regolith-linux/regolith-ftue/blob/2b2c1a9d73e93ebc67640876839ed853130f0e8e/regolith-init-term-profile) script also needs to use `trawlcat` as resource getter. This script is responsible for changing the terminal profile theme.

### Regolith-inputd

Similar to outut management, input management in wayland is also dependent on the compositor. This is the reason why changing input settings from `gnome-control-center` doesn't work out of the box. In order to get the input panel of gnome-control-center to work, I created [`regolith-inputd`](https://github.com/regolith-linux/regolith-inputd), which monitors the `gsettings` schema [`org.gnome.desktop.peripherals`](https://github.com/regolith-linux/regolith-inputd/blob/78349e1c9a77e4edcffbe5c245f8c58e138f90ea/org.gnome.desktop.peripherals.gschema.xml). This allows to map changes in input properties to sway configuration commands which are then applied using the `sway` IPC interface.

## What's next?

You can track work related sway integration in [this](https://github.com/orgs/regolith-linux/projects/8) GitHub projects page. All contributions are more than welcome. I have really enjoyed working on this project and I plan to continue even after I'm done with the final evaluation for Google Summer of Code.

## What doesn't work

- Lock Screen (Alternatives - swaylock)
- GNOME OSD overlays (Alternatives - wob)
- GNOME flashback / gnome-shell-related stuff
- Using any Display Manager other than GDM 3 causes issues
- swaybar status (requires several reloads to work properly)
- ilia keybindings page
- `gnome-session-quit` popup
