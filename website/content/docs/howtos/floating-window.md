---
title: "Launch an app in a Floating Window"
linkTitle: "Launch an app in a Floating Window"
description: >
  Most GUI apps can be launched in floating window mode.
---

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

Most X11 applications support a flag `--class` to specify the `class` under which the application runs.  Regolith is configured such that apps with a class of `floating_window` will launch with i3's floating window mode.

## Examples

### `gnome-terminal`

Launch the terminal in floating mode:

```bash
gnome-terminal --class=floating_window
```

### Firefox

Launch Firefox in floating mode:

```bash
firefox --class=floating_window
```
