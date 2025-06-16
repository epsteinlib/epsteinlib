---
title: "Verify behavior against upstream"
linkTitle: "Verify behavior against upstream"
description: >
  This page describes how to root cause the origin of bugs.
next: /docs/reference
---

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

Regolith is composed of several existing open-source projects plus some additional custom packages.  Sometimes when a bug is found in Regolith, actually it's a bug in an upstream project.  Determining the source of bugs is a first step in getting them fixed.

## Check `gnome-flashback`

Most of Regolith's GNOME support comes from the `gnome-flashback` project.  To run a `gnome-flashback` session, install the session and log back in to your account after selecting that session from the login screen.

1. Install the `gnome-session-flashback` session:

```bash
sudo apt install gnome-session-flashback
```

2. Restart your computer and then select the `GNOME Flashback (Metacity)` session as shown here:

{{< img "images/gnome-flashback-session-login.png" "">}}

3. Verify your behavior against Regolith.  If the problematic behavior also exists in `gnome-flashback` then the bug is likely there.  If the bug does not exist in the `gnome-flashback` then it's likely to be in Regolith.
