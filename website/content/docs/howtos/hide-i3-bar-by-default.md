---
title: "How to Hide i3-bar by default"
linkTitle: "How to Hide i3-bar by default"
description: >
  Hide i3-bar by default
---

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

Add the line
```yaml {filename="~/.config/regolith3/Xresources"}
wm.bar.mode: hide
```
to `$HOME/.config/regolith3/Xresources` file.
The next time you log in the bar should be hidden by default.
