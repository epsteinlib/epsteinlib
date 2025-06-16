---
title: "Disable Mouse Hiding (unclutter)"
linkTitle: "Disable Mouse Hiding (unclutter)"
description: >
  Disable unclutter-xfixes
---

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

By default the mouse cursor disappears if unused for a few seconds. If you prefer that the mouse cursor always remains visible on the screen, follow these steps:

1. Specify an empty value for the unclutter launcher program with an Xresources override:
  ```bash
  echo "i3-wm.program.unclutter: :" >> ~/.config/regolith3/Xresources
  ```
2. Log out and back in for changes to take effect.

