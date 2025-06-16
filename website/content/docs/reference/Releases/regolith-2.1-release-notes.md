---
title: "2.1 Release Notes"
linkTitle: "2.1 Release Notes"
weight: 89
description: >
  Release notes for Regolith 2.1
---

Regolith 2.1 is a minor release including several improvements, bug fixes, and optimizations.

## New Features

* Regolith Linux 2.1 ISOs (based from Ubuntu 22.04)
* Looks `dracula`, `ayu`, `ayu-dark`, `ayu-mirage`, `i3-default`, and `lascaille` added
* Looks can be loaded from user directories (`~/.config/regolith3/looks`)
* [`i3-swap-focus`](https://github.com/olivierlemoal/i3-swap-focus) is now available. Navigate to last window with keybinding `<super> .`

## Bug Fixes

* Ilia's file search dialog was broken on Jammy, this has been fixed
* `PtrSc` hotkeys working as in 1.6
* Numerous small fixes listed below

## Known Issues

* No localization for `ilia` 
* Inconsistent GTK theme in some dialogs in Jammy (`gnome-flashback` [fix in progress](https://bugs.launchpad.net/ubuntu/+source/gnome-flashback/+bug/1979120), will be available from Ubuntu repos)

## Bug Fix Releases

### 2.1.1

The Regolith Linux ISO was rebuilt as 2.1.1 due to a [packaging mistake](https://github.com/regolith-linux/regolith-desktop/issues/692).  The md5sum that may be used to verify the download was incorrect.  The 2.1.1 release fixes this by including the correct md5sum file.  While the 2.2.1 ISO was freshly built, it includes no package changes over 2.1.


## Changelog Delta from Regolith Desktop 2.0 to 2.1

See the [GitHub 2.1 release](https://github.com/regolith-linux/voulage/releases/tag/pkgbuild-20220718_013033) for additional details.

### Changes in `regolith-lightdm-config`:
```
b8d4503 Add attribution for NASA Juno image
fdbd150 Switch to Papirus icon theme. Set regolith as default session.
0910bd8 Add all configuration items to slick-greeter file
e95bd47 move from gtk greeter to slick. Update background
```

### Changes in `regolith-session`:
```
d1e21b9 Remove old look comment declarations
d2bee02 Add ability to stage looks in user directory ~/.config/regolith3/looks
ac6eb70 Check if user's directory does not exist before appending
```

### Changes in `i3-wm`:
```
937549f2 Merge tag '4.20.1'
```

### Changes in `i3-gaps-wm`:
```
b451b2a6 Merge tag '4.20.1'
```

### Changes in `regolith-i3-config`:
```
3f12389 Package changes
1e4f88c Switch to super-. for swap-focus based on user feedback
c5b43ff Package changes
6881ba4 Switch swap-focus to easier keybinding
6ab7598 Package changes
47c0545 Fix copypasta in i3-swap-focus
dff03ff Package changes
7f9ae0b Add new integration i3-swap-focus
647a658 Package changes
d7aeb82 Merge pull request #6 from jjsearle/default-border-style-option
df223f1 allow titlebar to be enabled via Xresources
```

### Changes in `regolith-look-default`:
```
3fba280 Reduce dependencies in default look, conform to adwaita
24ffab7 move to hicolor icon theme for default to reduce overhead
c187a2f Use colors in workspace labels
06841a8 Set titlebar style settings
589b264 Specify glyph for app laucher. Better glyph for help
```

### Changes in `regolith-look-extra`:
```
e0e049a Fix paths to compositor
a7018ec Decouple picom from looks. relates to https://github.com/regolith-linux/regolith-look-extra/pull/8
1d33fb1 Merge pull request #8 from mfisher87/consider-picom-deps
210c2d3 Merge pull request #6 from regolith-linux/nevil-icon-theme-default-switch
f70c4de Remove picom as a dependency of Nord and Nevil looks
632c175 updated Debian package metadata to install moka-icon-theme
1287be9 Add missing regolith-look-ayu-mirage look to debian packaging
13d125d Update looks to include lockscreen settings
17afcb7 Port lascaille from default look
7b5fb2e Add i3-default look
8ac4b27 Add window border styles to existing looks
883e434 Changed default nevil icon theme to Moka
493c746 Update workspace glyph icons across a few looks for more flair
b63b75a Add Ayu, Ayu-Mirage & Ayu-Dark
a600730 Specify glyph for app laucher. Better glyph for help
b3e6295 Add dracula look by @jjsearle
92f844a s/nord/nevil in composer-init
```

### Changes in `i3-swap-focus`:
```
c04ad13 Fixes for newer deb build tools
ccb5ce6 Add debian packaging
```

### Changes in `regolith-i3xrocks-config`:
```
11c7516 replace reboot command with saner default (#134)
086d0a2 Add control entry for i3xrocks-update
1d09bb8 Add tailscale IP (#132)
7f70e0a Add support for wireguard interface/tunnel (#118)
38da9e1 Add support to temp script for Tctl (#124)
0fe5a12 Remove extra space in output (#130)
e06b955 Add updates block (#131)
56ed590 move dynamic block window-name to the left-most position due to horizontal variance
ce6cebd Remove unneeded script for app launch action
86abc8c Load ilia over remontoire for help action
d400913 Adding a fix to not show anything on empty workspace (#129)
25fd18d Add support for bond interface in net-traffic (#128)
```

### Changes in `regolith-compositor-compton-glx`:
```
7980eb1 Allow config to be specified by Xres. Relates to https://github.com/regolith-linux/regolith-look-extra/pull/8
```

### Changes in `regolith-compositor-picom-glx`:
```
ee673f1 Allow config to be specified by Xres. Relates to https://github.com/regolith-linux/regolith-look-extra/pull/8
```

### Changes in `ilia`:
```
84687f2 Sync project versions
149d4e6 Filter on both summary and keybinding text in keybinding page to address https://github.com/regolith-linux/ilia/issues/34
eaf40b6 Revert to using x-terminal-emulator over gdk in support of https://github.com/regolith-linux/ilia/issues/31
156d478 Clean up page key event handling. Exit on system key events to address https://github.com/regolith-linux/ilia/issues/24
c0e4bd1 Fix double filtering in tracker page. Address https://github.com/regolith-linux/ilia/issues/23
3948bc8 Cleanup icon handing across pages. Fix icon handing for window view
52d350c Read floating node members of i3 tree response.  Addresses https://github.com/regolith-linux/ilia/issues/21
64943d8 Remove unused dependency, causing errors on ISO release
157df07 Remove default sorting of items from textlist page
```

### Changes in `regolith-default-settings`:
```
366893e Update default settings to match default look
```

### Changes in `remontoire`:
```
ab00f14 Fix version returned by program to match package version
```

### Changes in `plymouth-theme-regolith`:
```
b14a39f Use final logo image
```
