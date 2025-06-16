---
title: "3.2 Release Notes"
linkTitle: "3.2 Release Notes"
weight: 85
description: >
  Release notes for Regolith 3.2
---

## Regolith 3.2 Release Notes

Regolith 3.2 is a minor release focusing on Ubuntu 24.04 support, init system overhaul, and bug fixes.

## Features

* Ubuntu 24.04 (Noble) support
* `regolith-control-center` rebased from GNOME 46 (for supported releases)
* Sway session promoted from `beta` to `stable` status
* Overhaul of X11 and Wayland session initialization

## Deprecations

* Deprecating support for the following OS releases:
    * Ubuntu 20.04 (Focal)
    * Debian 11 (Bullseye)
* Deprecating support for Wayland session on Debian 12 (Bookworm) due to out-of-date build dependencies
* `regolith-control-center` has *not* been ported to GNOME 47 and so users on Ubuntu 24.10 will need to use `gnome-control-center` instead.


## Installation Instructions

For those using the `release-3_1` `apt` URL, update it to `release-3_2` and update your packages with `apt update` or similiar.  

For those using the "rolling" release label `release-current`, no change is necessary to get the update.

After `apt upgrade` is complete, verify that Regolith 3.2 is installed by:

```console
$ cat /etc/regolith/version 
REGOLITH_VERSION=3.2
```

Note: A reboot is required for the new version to load.

## Changes since Regolith 3.1

### Changes in `sway-audio-idle-inhibit`

```
2011f38 Merge branch 'upstream_main' into sync_upstream
a5b9689 fix: depend on sway instead of libwlroots11
```

### Changes in `ilia`:

```
b658923 chore: add layer_shell supported check
2e5029a fix: Ensure single instance using Gtk.Application
9c41eaf fix: deb tool config
353aff8 feat: launch apps and commands using systemd-run
```

### Changes in `picom`:

```
a959acc fix: update patch to enable config file
5c7f80e fix: update deb branch conf
9133602 fix: prevent picom autostart for regolith compat
df06cae fix: add missing dep
16731b1 fix: update dependencies for v11
d0ebbe8 chore: copy debian packging from v9 branch, add changelog entry
```

### Changes in `regolith-displayd`:

```
8c70872 chore: start daemon using systemd instead of regolith-displayd-init
f804063 fix: regolith-displayd first time startup failure via systemd
d72589e chore: let systemd manage regolith-displayd
```

### Changes in `regolith-lightdm-config`:

```
729dd65 fix: default session. feat: new wallpapper
```


### Changes in `regolith-look-extra`:

```
8ea5655 blackhole: Fix wallpaper.color.primary to use defined black
```

### Changes in `regolith-powerd`:

```
11684f8 chore: update Cargo.toml to match debian
db08cc4 feat: pause all audio players before suspend
c2bb7c3 chore: increase sleep time to account for loading of background image
2c1754c fix: also show the correct wallpaper when locking
b8f57f2 fix: use correct lock screen program when idling
```

### Changes in `regolith-sway-touchpad-gestures`:

```
e84d769 fix: require deprecated package to be removed before replacing
fcfabb1 chore: extract gesture partial to seperate source package for selective install based on target os due to compat issues on jammy
```

### Changes in `regolith-wm-config`:

```
1bdd1c4 fix: pull window binding definition out of partial, into root config, to prevent decl ordering issues
c5921b8 fix: specify valid color default. Load best quality image instead of worst.  Only set wallpaper if resolution match found. Cleanup.
44a562c feat: support overriding more keys via resources
66d2f1f fix: specify defaults for trawl. exit on error.  fix errors
29a4c50 fix: import I3SOCK into the systemd environment for sway
8511c1c chore: improve error handling
59cc3ad fix: use gtklock from path
2842af1 fix: update pkg name for touchpad package. remove bin package decl. add soft dep on grimshot.
1f45a52 fix(clamshell): directly execute actions
770bfe7 fix: pull nav binding definition out of partial, into root config, to prevent decl ordering issues
25ff6f1 feat: add 3 bindings for window and screen grabs with grimshot
566714c fix: ensure lid action is applied as per current power state
1ab9c54 fix: ensure lid action is applied as per current power state
824c5f9 feat(sway): add touchpad gestures support
b1f83c7 feat: use systemd-run for terminal and browser
30955c4 feat: customise lid close behaviour based on system state
3cf993c fix: add essential packages to recommends
```

### Changes in `trawl`:

```
6503886 feat: use absolute versioning for dependencies
```

### Changes in `regolith-wm-config`

```git
56a6e8d fix: use systemd to init xdg-desktop-portal
ef17348 fix: import dbus activation ennvironment before other config partials
24e8795 fix: remove redunant call to dbus-update-activation-environment
a2832ce fix: specify 'no-startup-id' to avoid waiting mouse cursor in X11, unclear of effect in sway, adding for consistency.
08cc043 feat: notify systemd about successful session activation
6d1ba64 Fix Mediakeys binding to support muting the mic
e81a80a feat: add package for running mate polkit under xwayland
a0e3c3c fix: mate polkit binary location is distribution dependant
```

### Changes in `regolith-desktop`

```git
5ce67ee chore: version bump to 3.2
```

### Changes in `i3status-rs`

```git
db13009 fix: disable battery to avoid global error on bar on targets without batteries on Ubuntu 22.04
```

### Changes in `xdg-desktop-portal-regolith`

```git
08cd870 fix: explicitly assign problematic portal interfaces (that fallback to GNOME) to regolith portal
```

### Changes in `regolith-displayd`

```git
b89ac56 amend version based on recent change
bd474a5 fix: use systemd to synchronize startup for displayd and kanshi (#5)
```

### Changes in `regolith-inputd`

```git
34cb76c feat: add systemd service
8bbbee3 fix: rename var delay -> sources
```

### Changes in `regolith-powerd`

```git
05afc0a feat: prevent logind actions for power key and lid
4f57160 feat: add systemd service
```

### Changes in `regolith-session`

```git
1f59e3e fix: add bookworm to legacy path
cc0020e fix: invert predicate on matching os release to gnome-session version
c79e188 feat: migrate sessions to systemd
2109934 fix: dont start session service starting targets
ab43221 chore: rename sway bin based on changes in sway-regolith package
```

### Changes in `regolith-wm-config`

```git
25c83a2 fix: add exec flag to fix spinner mouse cursor not present upon desktop load on x11
21774ac feat: restart regolith-init-kanshi service as user
9f68d28 feat: Update package dependencies version
148a376 feat: remove manual startup of gsd-integrations
3004ebe chore: manually start regolith-wayland.target if session doesnt use systemd
8b25827 feat: add keybindings to cycle between keyboard layouts
ba7841a chore: add packages for dbus-activation configs to Depends to ensure they're always installed with the session
5bc7cea fix: import XDG_CURRENT_DESKTOP into systemd_environment to ensure correct portal config is selected
```

### Changes in `sway-regolith`

```git
5ef388bf fix: remove desktop session file
b1857f95 fix: declare incompat w/ upstream
88b5271a fix: clean up package metadata
e88844b2 chore: update patches to use upstream bin name
6bcc07b7 chore: cleanup unrequired deps from the sway-regolith package (#22)
```
