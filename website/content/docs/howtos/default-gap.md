---
title: "Change the Default Window Gap"
linkTitle: "Change the Default Window Gap"
description: >
  Change the pixels between windows.
---

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

Like most configuration in Regolith, the i3-gaps gap size can be tuned via `Xresources` ([See here]({{< ref "/docs/Reference/xresources.md" >}}) for all Xresource definitions).

1. Create or add the following value to your `~/.config/regolith3/Xresources` file:

```yaml {filename="~/.config/regolith3/Xresources"}
wm.gaps.inner.size: 20
```

In this example we're setting the default gap to 20 pixels.

2. Reload the `Xresources` configuration:

```bash
regolith-look refresh
```

## Further Reading

See the [reference page for configurations]({{< ref "/docs/Reference/configurations.md" >}}) for more details about config files in Regolith.
