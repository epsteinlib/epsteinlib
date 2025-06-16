---
title: "Change the Terminal Font Size"
linkTitle: "Change the Terminal Font Size"
description: >
  Change the font used in the terminal
---

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

1. Create or add the following value to your `~/.config/regolith3/Xresources` file:

```yaml {filename="~/.config/regolith3/Xresources"}
gnome.terminal.font:	JetBrains Mono 14
```

2. Reload the `Xresources` configuration:

```bash
regolith-look refresh
```

3. Open a new terminal to see the change take effect.
