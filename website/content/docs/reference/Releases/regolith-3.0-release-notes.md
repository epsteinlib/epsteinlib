---
title: "3.0 Release Notes"
linkTitle: "3.0 Release Notes"
weight: 87
description: >
  Release notes for Regolith 3.0
---

## Regolith 3.0 Release Notes

[Regolith 3.0](https://github.com/orgs/regolith-linux/projects/26) is a major release focusing on a new [Wayland](https://github.com/orgs/regolith-linux/projects/8)-based session, Debian 12, and Ubuntu 23.04 support.

## New Major Features

* Ubuntu 23.04 (lunar) support
* Debian 12 (Bookworm) support
* Alpha release of Wayland session based on the [Sway compositor](https://swaywm.org/) (See [migration guide](https://regolith-desktop.com/docs/howtos/install-sway/))

## Migration Guide

{{< callout type="warning" >}}
Important information for existing Regolith users
{{< /callout >}}

1. The directory for user-staged configuration files has changed from `~/.config/regolith2` to `~/.config/regolith3`.
2. To generalize across X11 and Wayland sessions, Xresource keys have been renamed to remove "i3-".  This means that users with XResource overrides on earlier versions of Regolith will need to update the key names.  Specifically, keys starting with `i3-wm` now start with `wm`.  For example, in Regolith 2.x to change gaps size, `i3-wm.gaps.inner.size` is being renamed to `wm.gaps.inner.size` in Regolith 3.0.
3. Some functionality may now require you to install certain packages manually. An example here is the `regolith-wm-user-programs` package that runs programs defined in your Regolith configuration declared by `wm.program.<N>:`.

## Changelog Delta from `beta 3` to `beta 4`

### Changes in `regolith-displayd`:
```
4e2606b feat: add fractional scaling for resolutions greater than 1920 x 1080
0ea712d chore: linting changes
```

### Changes in `regolith-i3xrocks-config`:
```
3eed9da Package changes
7f25895 Use a … instead of ... to show more of the artist and title (#143)
2f73725 Use value colors for battery percentage and charge (#148)
```

### Changes in `regolith-session`:
```
37f0773 fix: add conflicts to avoid collision on ext package
c0bf9b8 move autostart file to new package
76d8d1d fix: move idle autostart file to seperate package as soft dep to maintain compat w/ upstream flashback
c072b4f fix: path to idle monitor daemon incorrect
e092d1d fix: Ship idle monitor autostart desktop file.  Addresses https://github.com/regolith-linux/regolith-session/pull/15
7f76fdf fix: install picom compositor by default on x11 session
59120ee feat: enable loading of last installed theme on first session install in xorg session
fe33317 fix: changelog version invalid
b8c6564 fix: on first lanch of sway session, call load_look() function to setup environment
9c9c42b add error handling in case that Xres/trawl look path is invalid
f67e8ac feat: if a look is not set in xres/trawl, load the most recently installed look instead of default
0389967 fix: push regolith-default-settings down to x11 session, possibly unneeded for wayland (experiemental)
```

### Changes in `regolith-wm-config`:
```
f5892dc fix: inconsistent behavious of XF86AudioPlay
ba26d1f chore: update homepage
d85be5a fix: add i3 gaps to debian control
```

### Changes in `regolith-desktop`:
```
821c31c fix: push regolith-default-settings down to x11 session, possibly unneeded for wayland (experiemental)
```

## Changelog Delta from `beta 2` to `beta 3`

### Changes in `regolith-wm-config`:
```
d95e330 fix: inconsistent behavious of XF86AudioPlay
9d98515 fix: add i3 gaps to debian control
32fd6d4 Remove deprecated dependency i3-gaps
21da9b6 Fix conflict between deprecated regolith-wm-ilia and wm specific updates
6e697ad fix: add missing package declaration in control
05b8e9d fix: use 'bindsym --to-code' on sway for more consistent key mapping behavior with i3
```

### Changes in `regolith-session`:
```
e0012a4 feat: add additional diagnostic infor about the current session
```

### Changes in `ilia`:
```
62aef17 fix(wayland): ilia returns empty keybindings page if a directory included in the configs doesn't exist
c3a03f8 more window manager agnostic renames
8cfc7e3 rename I3Client to IPCClient
c78fa1b rename read_i3_config to read_config
```

### Changes in `regolith-i3xrocks-config`:
```
7f25895 Use a … instead of ... to show more of the artist and title
2f73725 Use value colors for battery percentage and charge
```

## Changelog Delta from `beta 1` to `beta 2`

### Changes in `regolith-wm-config`:
```
445edee feat: add common directory for shared configs between i3 and sway
d3c5ca2 fix(#866): incorrect path for i3 user config directory
```

## Changelog Delta from Regolith Desktop 2.2 to 3.0

### Changes in `regolith-session`:
```
41dbd95 Making regolith-look-default hard dependency
f17a946 Migrate to wm agnostic packages in x11
c5e7900 version bump
f5a9901 Minor cleanup of package metadata
14e41f7 Changed source package format
3f042ff Add ilia and rofication and associated config to sway session recommends for usable setup
765ba57 Add dependency on virt look package to ensure resources are loaded
a98da23 Change package format
385c478 Fixed workspace naming for resource cleanup
e78230a Added gnome-terminal to recommends
c4e17b3 Prevent loading the incorrect look
49b09f7 Added gnome-settings-daemon as dependency for regolith-session-sway
9968419 Move desktop dep to flashback session to keep desktop wm indenpendent
5aba956 cleanup
5f9f1d0 Update package deps for wayland
32f9eb9 Use quilt packaging format in hopes of resolving build error
be04b61 Update package metadata, attempt fix of regolith-diagnostic confict
5e80731 Fix version so package will build into repo
a6b7daa Send tick before reload on sway
c1bdf0f Remove RfKill from required components
4c5aa9d Add regolith-sway-config as dependency
71d906a Updated depends for regolith-session-sway
9ee4d88 Fix look not loading on startup
03ab9a4 Load environment variables from systemd-environment-d
933ccbf Fixed typo
9c9a33a Added support for unsupported GPU
```

### Changes in `i3-wm`:
```
b85da284 release i3 4.22
```

### Changes in `regolith-look-default`:
```
d1f8351 Update upstream theme package for lunar
403c89a Remove i3 from wm filename
57af57f Refactor xres namespace to be x11 agnostic.  i3-wm -> wm
```

### Changes in `regolith-look-extra`:
```
e7c9247 Remove i3 from wm filename
8076d02 Refactor xres namespace to be x11 agnostic.  i3-wm -> wm
ad3ec5f Lighten borders of inactive/unfocused tabs
d33d503 Add README placeholders for testing workflow
7d78827 Style unfocused/inactive tabs the same way as inactive workspaces
0f88ed5 Revert change to original theme palette
67c946e Use darker color for nord base02
cbfd4e6 Cleanup deprecated picom options
a6ea33d Misc tweaks to extra looks
```

### Changes in `i3-next-workspace`:
```
fc40c1f Relax search on workspace name to work with Regolith 3.0
```

### Changes in `regolith-i3xrocks-config`:
```
e70a24b fix(keyboard-layout): open keyboard menu instead (#145)
```

### Changes in `regolith-desktop`:
```
e44a757 Update dependencies for support of x11 and wayland (#812)
291c3be Version bump for next release
```

### Changes in `regolith-ftue`:
```
017a50b Update sig for build
```

### Changes in `regolith-compositor-picom-glx`:
```
0d7ad22 Remove deprecated picom options
```

### Changes in `ilia`:
```
86f5694 Bump version to sync w/ jammy branch
b8ea9d6 add keybings to help
5e096d6 add missing vim keybinding to navegation
fe91a57 Fix all page mode curor navigation.
490ef9e Fix inconsistent widget focus init with new page function.
bd977ec Remove async app loading for https://github.com/regolith-linux/ilia/issues/51
```

### Changes in `regolith-rofication`:
```
399b528 Remove wm dependencies for wayland support
```

### Changes in `regolith-wm-config`:
```
cbe6f87 fix: incorrect arguments to hide_cursor
7b3149b fix: add missing control entry for regolith-sway-unclutter
7ce2e40 fix: add missing .install file for regolith-sway-unclutter
cf333da fix: XF86AudioPrev didnt play the previous track
e36deff Remove duplicate "Shift"
59d2dba Add regolith-sway-polkit to depends
ed5ccfa Packaging changes
88cca6b Fix gtklock binding not working when wallpaper not set
7d7550c Added regolith-sway-gtklock to recommends of regolith-sway-root-config
702d92d Remove unknown information fields
4ed1662 Added config for regolith-sway-audio-idle-inhibit
4f8c2e5 Fix typo
e42d548 Add expected packages to recommends
3a11528 Removed redundant variable initialization
c943500 Add locked flag to some keybindings
447498c removed self-conflict and self-replace for regolith-wm-control-center
e437d32 Added packages to regolith-sway-root-config 'Recommends'
8ca9ab3 Updated package version
11a12bf Replace source package regolith-i3-config with regolith-wm-config
f258e51 Remove regolith-wm-session dependency from regolith-wm-config binary package
fc23d0c Fixed workspace names on initial load, sway
9117eeb Config fixes
5b9b7fe Fixed packaging
daebaf4 Updated sway/config.d/81_gsd_integrations and added regolith-powerd as dependency
6b1e430 Exec regolith-powerd on wm_startup, regolith-displayd
324da56 Added config partial for gtklock
01127a5 Add config partials and i3status-rust config
46d09bb Fixed dependencies
9bf981f More whitespace cleanup in debian control
b81cde0 Debian package fixups.  Add changelog entry
c02c847 Renamed conf.d directories to config.d
f31cd23 Added packages to configure sway session
1905f46 Changed variable naming for sway configs
ccf73eb Add root config files and packages for the root config files for i3 and sway.
01b139c Update file packages based on feedback from PR.
d89da22 Initial import of debian and wm config from legacy regolith-i3-config repo
09fd304 Initial commit from regolith-i3-config
```

