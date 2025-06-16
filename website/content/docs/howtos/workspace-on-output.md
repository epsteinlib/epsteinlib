---
title: "Assign a workspace to a specific Output (Display)"
linkTitle: "Assign a workspace to a specific Output (Display)"
weight: 1
description: >
  Specify that a given workspace should always be displayed on a specific monitor
prev: /docs/howtos
---

## i3

{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

In order to configure i3 to assign a workspace to a given monitor, first find the device name for the monitor you wish to use with `xrandr`:

```console
$ xrandr --listmonitors
Monitors: 1
 0: +*eDP-1 1920/309x1080/173+0+0  eDP-1
```

In this case, the output name is `eDP-1`. Next add a line like this into [a staged copy of your i3 config file]({{< ref "stage-configs">}}):

```text
workspace "$ws1" output DP1
```

Log back in for the changes to take effect.

## Sway

The concept is the same, but `swaymsg` is used to display the available outputs.

```console
$ swaymsg -t get_outputs
Output HDMI-A-2 'ViewSonic Corporation VG2448 V5E192121969' 
[...]
```

To assign workspace *$ws1* to the HDMI-A-2 output put the following line into [a staged copy of your sway config file]({{< ref "stage-configs">}}):

```text
workspace $ws1 output HDMI-A-2
```
