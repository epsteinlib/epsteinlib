---
title: "1.6 Release Notes"
linkTitle: "1.6 Release Notes"
weight: 91
description: >
  Release notes for Regolith 1.6
---

Regolith R1.6 is a feature release which includes several improvements and optimizations. Regolith 1.6's primary feature is a return of the ISO installation method.

## Features

| Feature                                 | Screenshot                                                                                                                                                                                                           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Return of Regolith Linux ISO      | {{< img "images/releases/regolith-1.6.0-login.png" "Regolith 1.6 Login Screen" >}}    | Up to Regolith 1.4.1, we hand-built each ISO for major version releases. This became too brittle and time intensive, and so work began on [build automation to generate the ISO](https://github.com/regolith-linux/regolith-ubuntu-iso-builder).  This took quite awhile to complete but is [now available](https://github.com/regolith-linux/regolith-ubuntu-iso-builder/releases) for the 1.6 release.                                                                                          |                  |
| New Look: Midnight                     | {{< img "images/releases/regolith-look-midnight.png" "Midnight Look" >}} | Fabian (naib864) has provided a distinctive new look with new font, icons, and GTK theme. |
| New Look: Solarized Light         | {{< img "images/releases/regolith-look-solarized-light.png" "Solarized Light Look" >}} | Ricardo Pérez (ricpelo) has provided a solarized light theme.
| New Status Indicator: app-launcher| | Regolith user Laurent provided the idea to add a status indicator to launch Rofi, and so `i3xrocks-app-launcher` was added for those that prefer to use the mouse to launch apps. | 
| Upgrades to i3-gaps and Rofi | | [i3-gaps is updated to 4.19.1](https://github.com/Airblader/i3/releases/tag/4.19.1) and [Rofi to 1.6.1](https://github.com/davatorium/rofi/releases/tag/1.6.1). | 

## Fixes

Have a look at [the R1.6 project page](https://github.com/orgs/regolith-linux/projects/12) for a list of bug fixes.

## Changelog Delta from Regolith 1.5.3 to Regolith 1.6

```
dracula-gtk (1.1.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add icons to dracula theme

regolith-desktop (2.122-1bionic) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add app launcher to complete desktop metapackage. Fix unneeded dep on mutter-common.

regolith-desktop (2.118-1bionic) bionic; urgency=medium

  [ Ken Gilmer ]
  * Version bump to 1.6.0

regolith-ftue (1.1.2-1) bionic; urgency=medium

  [ Jordan Callicoat ]
  * More Xresources options for GNOME Terminal

regolith-i3-gaps-config (2.8.6-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * specify default workspace layout mode

regolith-i3-gaps-config (2.8.5-1) bionic; urgency=medium

  [ Bence Ferdinandy ]
  * adding $mod+s for sound setting (#36)

  [ Ken Gilmer ]
  * Fix urls to howtos. Relates to https://github.com/regolith-linux/regolith-i3-config/pull/1.

regolith-i3xrocks-config (3.6.4-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove invalid Xresource key segment from microphone block

regolith-i3xrocks-config (3.6.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add missing dependency to net-traffic script on iproute2

regolith-i3xrocks-config (3.6.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix control file
  * Fix config for app-launcher, remove incorrect interval on next-workspace

regolith-i3xrocks-config (3.6.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Take nicolas-lopez's contribution via gist as microphone status indicator.  Related to https://github.com/regolith-linux/regolith-i3xrocks-config/issues/105.
  * Add app launcher status indicator

regolith-rofication (1.3.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove debug logs that show up in /var/log/syslog (#9)

regolith-rofication (1.3.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Require i3xrocks-rofication to be upgraded first, due to file collision between i3xrocks-rofication and regolith-rofication

  [ Carles Sala ]
  * Close window after deleting the last notification (#8)

regolith-styles (2.9.7-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Rename midnight look background file to avoid collision w/ ubuntu-wallpapers-focal.

regolith-styles (2.9.6-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Another attempt to fix pop-os wallpaper reference

regolith-styles (2.9.5-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix pop os gnome settings

regolith-styles (2.9.4-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix pop-os background. Change image filename to avoid conflicts

regolith-styles (2.9.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add background image missing from midnight theme

regolith-styles (2.9.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix color value unparsable by i3

regolith-styles (2.9.1-1) bionic; urgency=medium

  [ Ricardo Pérez ]
  * Fixes Gtk theme (SolArc-Dark -> SolArc) (#44)

regolith-styles (2.9.0-1) bionic; urgency=medium

  [ Ricardo Pérez ]
  * New solarized-light theme (#43)

  [ Ken Gilmer ]
  * Package solarlized-light look contribution from @ricpelo
  * Add icon for app launcher shortcut

regolith-styles (2.8.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add app launcher glyph

regolith-styles (2.8.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix font specification for midnight look

regolith-styles (2.8.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * midnight theme package fixes

regolith-styles (2.8.0-1) bionic; urgency=medium

  [ Fabian ]
  * Add midnight theme (#37)

  [ Ken Gilmer ]
  * Add dependencies to midnight look

regolith-styles (2.7.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Specify different color for indicator to make window orientation indicator visible in dracula.  Fixes https://github.com/regolith-linux/regolith-styles/issues/40.
  * Revert to single column mode in lascaille
  * Use accent color for orientation border decoration in lascaille

regolith-styles (2.7.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Rofi app launcher in two columns

regolith-styles (2.6.17-1) bionic; urgency=medium

  [ HieuTranTH ]
  * dracula: Fix a typo in i3-wm config file (#41)
  * ubuntu: Fix color_base2 incomplete hex code (#42)

  [ Ken Gilmer ]
  * Add highlight to bar for non-active focused windows for better multi-monitor experience (lascaille)

regolith-system (1.5.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add recommends for Regolith system packages

regolith-system (1.5.2-1) hirsute; urgency=medium

  [ Ken Gilmer ]
  * remove deprecated dependency

rofi (1.6.1-1regolith1) bionic; urgency=medium

  * Adapted package metadata from Ubuntu for 1.6.1 release 

ubiquity-slideshow-regolith (168.5-regolith1) focal; urgency=medium

  [ Ken Gilmer ]
  * Update virtual package dependency for focal

```
