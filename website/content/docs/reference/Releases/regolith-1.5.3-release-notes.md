---
title: "1.5.3 Release Notes"
linkTitle: "1.5.3 Release Notes"
weight: 92
description: >
  Release notes for Regolith 1.5.3
---

Regolith 1.5.3 is a bug fix release.  The following list contains brief descriptions of fixes and updates per package since the Regolith 1.5.2 release:

## Changlogs

### Release Notes for i3xrocks

i3xrocks (1.3.5-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Do not print system file error if user files do not exist.  Addresses https://github.com/regolith-linux/regolith-i3xrocks-config/issues/114.


### Release Notes for regolith-desktop

regolith-desktop (2.111-1bionic) bionic; urgency=medium

  [ Ken Gilmer ]
  * Move from i3xrocks-todo to regolith-todo to address https://github.com/regolith-linux/regolith-todo/issues/6.  Version bump in version file for 1.5.3.

### Release Notes for regolith-i3xrocks-config

regolith-i3xrocks-config (3.5.7-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * remove debugging output flag in net-traffic script

### Release Notes for regolith-styles

regolith-styles (2.6.16-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Tweak to width to better match pixel equivalent at 1080p.


regolith-styles (2.6.15-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Use relative width over pixels for better HiDPI support of Rofi dialog.  Idea from @eduncan911 in https://github.com/regolith-linux/regolith-desktop/issues/417.

### Release Notes for regolith-todo

regolith-todo (1.2.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Specify that regolith-todo replaces i3xrocks-todo.  Addresses https://github.com/regolith-linux/regolith-todo/issues/6.


### Release Notes for remontoire

remontoire (1.4.1-1) bionic; urgency=medium

  [ Mark Stosberg ]
  * Document `-i`, `-p`. (#11)

  [ Michael Vetter ]
  * Add openSUSE install instructions

  [ Ken Gilmer ]
  * Add fonts-font-awesome as soft dependency based on reference in CSS file and user feedback.  Addresses https://github.com/regolith-linux/remontoire/issues/16.
  * Cleanup


