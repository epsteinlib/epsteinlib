---
title: "Customize Compositor"
linkTitle: "How to Customize Compositor"
description: >
  How to customize Regolith's compositor
---

{{< callout type="info" >}}
Pay special attention to this section if you are experiencing odd visual
glitches or slow graphics performance.
{{< /callout >}}

A compositor is a UI component that applies visual effects to windows before
they are rendered on-screen. Many desktop environments integrate a compositor
directly into the window manager, making it difficult to switch out or disable.
In Regolith, the compositor is defined as a pluggable “extension point” in the
packaging system. This means that compositors can be switched out simply by
installing the packages that contain them. The underlying packaging system will
ensure there are no conflicts and that all the dependencies of a given
compositor are also installed.

{{< callout type="warning" >}}
In Regolith 2.x, the default compositor is "no compositor". In
Regolith 1.x, the default was Picom.
{{< /callout >}}

## Finding available compositors

The following command will list all compositors configured to work with
Regolith:

```bash
apt search regolith-compositor-
```

You should find at least the following:

- `regolith-compositor-none` **[default]**: No compositor. Best performance, no
  visual effects.
- `regolith-compositor-picom-glx`: Uses `picom`. Recommended for most users to
  resolve screen tearing and add visual effects. Modern and well-maintained
  fork of `compton`. Was the default compositor in Regolith 1.6.
- `regolith-compositor-compton-glx`: An ancestor of `picom` which is no longer
  maintained.
- `regolith-compositor-xcompmgr`: An early compositor that may offer more
  consistent performance if more complex compositors are problematic.

## Installing a compositor

{{< callout type="info" >}}
Due to the way that the compositor is managed by i3-wm, **you must log out and
back in for compositor changes to take effect**. Restarting i3 is not
sufficient.
{{< /callout >}}

Run the following package install command to swap out the compositor. In this
case we’ll install `picom`:

```bash
apt install regolith-compositor-picom
```

**Log out and then log back in** and the new compositor should be loaded. You
can verify by checking the process list for `picom`:

```bash
ps aux | grep picom
```

## Overriding Regolith default compositor configuration

### Compton/Picom

To provide your own Compton/Picom compositor config, copy the default or create
your own and save it as `~/.config/regolith3/picom/config`. Upon next session it
will be loaded instead of the default config `/etc/regolith/picom/config`.
