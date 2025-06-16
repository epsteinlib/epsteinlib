---
title: "Change the lock screen program"
linkTitle: "Change the lock screen program"
description: >
  Configure Regolith to use a lock screen other than `gnome-flashback`
---

## Configure `gnome-flashback` to Disable Built-in Lockscreen

By default `gnome-flashback` handles the lock screen for Regolith. This can be adjusted by modifying it's configuration via the `gsettings` command line program:

```bash
gsettings set org.gnome.gnome-flashback screensaver false
```

## Cause Your Lock screen program to Run from i3 Keybinding

Specify the following override line in your `~/.config/regolith3/Xresources` file (substituting `your-script.sh` with the command you wish to use for the locking your screen):

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

```yaml {filename="~/.config/regolith3/Xresources"}
wm.program.lock: your-script.sh
```

For example, assuming you have `i3lock` installed:

```bash
echo "wm.program.lock: /usr/bin/i3lock" >> ~/.config/regolith3/Xresources
```
