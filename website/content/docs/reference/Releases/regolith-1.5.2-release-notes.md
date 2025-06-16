---
title: "1.5.2 Release Notes"
linkTitle: "1.5.2 Release Notes"
weight: 93
description: >
  Release notes for Regolith 1.5.2.
---

Regolith 1.5.2 is a bug fix release.  The following list contains brief descriptions of fixes and updates per package since the Regolith 1.5.1 release:

## Changlogs

### Release Notes for regolith-compositor-picom-glx

  * Remove manual fade flag passed to compositor in favor of specification in config file.  Allows users to specify fade parameters via config file.

### Release Notes for regolith-ftue

  * Check update flag before changing gnome-terminal.  Allows users to turn gnome-terminal profile generation off.


### Release Notes for regolith-i3xrocks-config

  * Added gawk specific testing for temp script (#108)
  * Filter all interfaces associated with VPN traffic while determining default route (#110)
  * Capture error cases involving invalid network device paths and sleep before returning to resolve CPU consuming loop due to  interval type of i3blocks.
