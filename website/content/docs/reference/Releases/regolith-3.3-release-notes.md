---
title: "3.3 Release Notes"
linkTitle: "3.3 Release Notes"
weight: 84
description: >
  Release notes for Regolith 3.3 (Beta 1)
prev: /docs/reference/releases
---

## Regolith 3.3 Release Notes (Beta 1)

Regolith 3.3 is a minor release focusing on Ubuntu 25.04 support, numerous build system improvements, and sway 1.10.

## Features

* Ubuntu 25.04 (Plucky) support
* Ubuntu 24.10 (Oracular) support (without `regolith-control-center`)

## Known Issues

* `regolith-control-center` has *not* been ported to GNOME 47 (but has for GNOME 48) and so users on Ubuntu 24.10 need to use `gnome-control-center` for all system configuration tasks.

## Installation Instructions

### Beta 1

The `apt` URL for 3.3 beta 1 for Ubuntu 25.04 (Plucky) is:

```
deb [arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] https://archive.regolith-desktop.com/ubuntu/testing plucky main
```

Note the `amd64` architecture.  For ARM-based systems, change that to `arm64`.  For Ubuntu 24.04 (Noble) it is:

```
deb [arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] https://archive.regolith-desktop.com/ubuntu/testing noble main
```

Refer to the [installation instructions](/docs/using-regolith/install/) to update the line for other Debian variants.

Installation on a freshly installed Ubuntu or Debian system can follow the standard installation instructions.  When upgrading from a previous release
using `do-release-upgrade`, the Regolith package repo has to be re-added to your system as the installer removes it.  Ensure the following packages are installed after the Ubuntu release process has completed and the Regolith package repository has been added back to your `apt` configuration:

* `regolith-i3-ilia` and/or `regolith-sway-ilia`
* `regolith-i3-control-center-gnome` and/or `regolith-sway-control-center-gnome` on Ubuntu 24.10 only.

After `apt upgrade` is complete, verify that Regolith 3.3 is installed by:

```console
$ cat /etc/regolith/version 
REGOLITH_VERSION=3.3
```

## Changes since Regolith 3.2

### Changes in `i3-next-workspace`

```
b38150d fix: workspace regex
```

### Changes in `ilia`

```
0c84c06 Merge pull request #105 from manthanabc/clipboard
6bc2cc7 fix: modified ctrl-v past behavior
be6542c enhanced searchBox of ilia
e8d8678 updated the deprecated meson.source_root
b49f0f8 changed to not use nohup, but execute command directly (output is shown by invoking shell)
ceb3dbe changed command execution from bash to nohup (part of coreutils)
66a7c42 added -q option to --help message
c9bd471 added -q option to CommandPage.vala to allow quiet execution (i.e., no x-terminal-emulator is started and output is discarded)
```

### Changes in `picom`

```
4235b67 feat: use newer version of libpcre which is 2 not 3!
```

### Changes in `regolith-archive-keyring`

```
b57b588 feat: d/ initial commit
b6cac12 feat: add archive-key.asc file
9b0f177 Initial commit
```

### Changes in `regolith-control-center`

* New port from GNOME-48

### Changes in `regolith-desktop`

```
200b6fc chore: prepare for regolith v3.3 release (#1106)
9a08e39 chore: remove trailing empty line from version file
```

### Changes in `regolith-session`

```
d2709e2 fix: remove deprecated package from control file. gnome-flashback no longer requires it, has it's own policykit impl (#45)
38b14dd Fix/app centre (#41)
```

### Changes in `regolith-wm-config`

```
da9f594 fix: cause conflicts on control center packages so only one can be installed at a time
7b44c91 chore: remove startup for gnome-flashback-media-keys
dbe5776 fix: add  gnome-control-center configs for sway
```

### Changes in `remontoire`

```
7511383 Update main.vala
5ea9cd8 Fixing version on main.vala
```

### Changes in `sway-regolith`

* Ported Sway v1.10 to `sway-regolith`

### Changes in `trawl`

```
c879c18 fix: update build files for changed service file location
051214c Fix typos in README.md
```

### Changes in `ubiquity-slideshow-regolith`

```
5381645 feat: overall redesign of regoilth slides
84c7f94 chore: upgrade to version 202 from upstream
```

### Changes in `xdg-desktop-portal-regolith`

```
cd012f8 fix: build regression on plucky.  Seems like meson cannot find correct dependency without updating debian package dependency to slightly different name (#13)
3d73f37 feat: Implement set_wallpaper_uri method in Wallpaper.Portal (#3)
67214d1 Fix/jammy build (#6)
```
