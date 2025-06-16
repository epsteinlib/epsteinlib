---
title: "3.1 Release Notes"
linkTitle: "3.1 Release Notes"
weight: 86
description: >
  Release notes for Regolith 3.1
---

## Regolith 3.1 Release Notes

[Regolith 3.1](https://github.com/orgs/regolith-linux/projects/26) is a minor release focusing on Ubuntu 23.10 support and bug fixes.

## Features

* Ubuntu 23.10 (Mantic) support
* Debian Testing support for Wayland session
* `i3-next-workspace` [ported to Sway](https://github.com/regolith-linux/childe)
* `regolith-control-center` rebased from GNOME 45 (for supported disto releases)
* Sway session promoted from `alpha` to `beta` status

## Installation Instructions

Follow the standard installation instructions, but substitute `testing` instead of `release-3_0`.  For example, this is the `apt` line for Ubuntu Mantic:

```bash
echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
https://regolith-desktop.org/testing-ubuntu-mantic-amd64 mantic main" | \
sudo tee /etc/apt/sources.list.d/regolith.list
```

After `apt upgrade` is complete, verify that Regolith 3.1 is installed by:

```console
$ cat /etc/regolith/version 
REGOLITH_VERSION=3.1
```

A reboot is required for the new version to load.

## Changelog Delta from Regolith 3.0

### Changes in `regolith-control-center`

```text
8545068e7 chore: add Regolith-Wayland in OnlyShowIn
08a8514db fix: disable whoopsie for debian compat
```

### Changes in `regolith-displayd`

```text
aa3d572 fix: cannot reenable monitors
c4ebc88 fix: disabled monitors should reflect in mointors instead of logical_monitors
5cdb6df feat: add support for builtin monitor
```

### Changes in `regolith-i3xrocks-config`

```text
320980c fix: relax dependency on control center.
```

### Changes in `regolith-session`

```text
dafa3b7 fix: add (soft) dependency to xdg-desktop-portal-regolith due to gnome 45 interop issues.
1a1c5e0 fix: set XDG_CURRENT_DESKTOP, overriding previous value if already set
8f5da99 fix: set XDG_CURRENT_DESKTOP for wayland session (fixes screencasting)
```

### Changes in `regolith-wm-config`

```text
ea53b18 fix: make all non-root config partials soft dependendencies for more flexibility
6bd920a chore: rename resources namespaced in regolith.clamshell.* to wm.clamshell.*
04eee9b feat: implement clamshell functionality for sway
db38632 fix: force removal of deprecated ilia integration package
7f1bae9 fix: enable carry mode in childe
b9d42e8 fix: prevent broken legacy package from remaining installed for next-workspace
6992986 fix: properly segment i3 and say versions of next-workspace partial
b36e821 feat: add next-workspace feature to sway session, copied from x11 i3-next-workspace
6356077 [15_base_launchers] add vm.program.browser
a777246 [15_base_launchers] add vm.program.terminal
9b1706d Move default-style to soft dep for both sessions.  demote 3 more for sway session to soft per PR feedback
d21c32c fix: remove (redundant) hard dependency from sway session to rofication. Soft dependency comes via regolith-wm-config.
```

### Changes in `sway-regolith`

```text
75311f9a update patch filename for consistency
b94bd17b fix: do not ship wallpaper, provided by sway-wallpaper package
6e56e53f fix: remove git versioning and fix c warning to appease debian builder on mantic
8dd63940 feat: upgrade packaging to sway v1.8.1
2a2e3b93 fix: include wl_protos_src in sway_sources
9160dab1 feat: apply patch for dbus tray menu
34f5363c Applied patch for regolith compatibility
227ac44a debian: import packaging from v1.7
```

### Changes in `sway-audio-idle-inhibit`

```text
162c01f chore: bump libwlroots dependency version
```

### Changes in `trawl`

```text
dfcabfb chore: pin zvariant version to 3.4.1
```

### Changes in `xdg-desktop-portal-regolith`

```text
* Initial implementation, forked from `xdg-desktop-portal-elementary`.
```

