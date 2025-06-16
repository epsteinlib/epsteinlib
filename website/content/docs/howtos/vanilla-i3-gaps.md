---
title: "Launch Default i3 Session"
linkTitle: "Launch Default i3 Session"
description: >
  Install package `i3-gaps-session` to keep using i3 without uninstalling Regolith
---

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

To do so on Ubuntu, run `sudo apt install i3-gaps-session` to add an option on login that starts vanilla i3-gaps without regolith.

## Use Case

This is useful for, say, running an application that does not integrate well with GNOME. For example, Spotify's unmaintained official Linux client mishandles its file selection dialogue and crashes when

1. Setting user profile photo
2. Setting album cover photo
3. Adding a local directory to the library paths

under GNOME and Regolith. Switching to vanilla i3 solves the crash.

_Source: [Ken's comment on the relevant issue on GitHub](https://github.com/regolith-linux/regolith-desktop/issues/21#issuecomment-729776535)._
