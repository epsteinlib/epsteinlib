---
title: "Support Multiple Batteries"
linkTitle: "Support Multiple Batteries"
description: >
  Display status of the batteries per instance
---

## i3
{{< callout type="error" >}}
This page was copied from the [Regolith 1.x website](https://regolith-linux.org) and has not been updated for Regolith 2.  It may contain out of date information.
{{< /callout >}}

Out of the box the primary battery status will be displayed in the i3-bar. For some laptops that have multiple batteries, e.g. ThinkPads, you want to have both batteries displayed:

1. [Stage your i3xrocks configuration files]({{< ref "add-remove-blocklets#customize-bar" >}}), then open `~/.config/regolith3/i3xrocks/conf.d/80_battery` in your editor of choice.
2. Find the line that pertains to the battery blocklet:

```toml {filename="~/.config/regolith3/i3xrocks/conf.d/80_battery"}
# Battery indicator
#
# Displays total percentage charge left of specified battery.
[battery]
interval=30
# instance=BAT0
```

3. List all the battery instances running on your machine by running:

```sh
ls /sys/class/power_supply/

# AC BAT0 BAT1
```

4. Update the staged config file to have a new blocklet as such:

```toml {filename="~/.config/regolith3/i3xrocks/conf.d/80_battery"}
...
# Battery indicator
#
# Displays total percentage charge left of specified battery.
[battery]
interval=30
instance=BAT0

[battery]
interval=30
instance=BAT1
```

## Sway

1. Stage your i3status-rust config file:
```bash
cp /etc/regolith/i3status-rust/config.toml ~/.config/regolith3/i3status-rust/config.toml
```

2. Configure to supply i3status-rust with your config file:
```bash
echo "wm.bar.status_config: ~/.config/regolith3/i3status-rust/config.toml" >> ~/.config/regolith3/Xresources
```

3. List all the battery instances running on your machine by running:

```bash
ls /sys/class/power_supply/

# AC BAT0 BAT1
```

4. Update the staged config file to contain blocks for both of the batteries:

```toml {filename="~/.config/regolith3/i3xrocks/conf.d/80_battery"}
[[block]]
block = "battery"
interval = 10
[...]
device = "BAT0"

[[block]]
block = "battery"
interval = 10
[...]
device = "BAT1"
```

## Further Reading

See the [reference page for configurations]({{< ref "/docs/Reference/configurations.md" >}}) for more details about config files in Regolith.
