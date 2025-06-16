---
title: "Install the Sway Compositor"
linkTitle: "Install the Sway Compositor"
description: >
  Upgrade to the Wayland-based Regolith session
---

To install the wayland session, follow these steps:

1. Install the `regolith-session-sway` package and it's dependencies.
2. Reboot
3. Select "Regolith on Wayland" session in the login panel

## Hardware Support

Those using proprietary nvidia drivers for their graphics cards must set `regolith.sway.unsupported_gpu: true` in their XResources. This is not an issue with the open source drivers nouveau.  See more [here](https://github.com/swaywm/sway/issues/5392).

## Bar Changes

The wayland session uses [`i3status-rs`](https://github.com/greshake/i3status-rust) rather than `i3xrocks`.  `i3status-rs` provides a crisper and more resource-efficient approach to system monitoring.  However at this time the `i3xrocks` packages for status indicators are not integrated with `i3status-rs`.  To use a custom configuration with `i3status-rs`, users will need to put their config file at the location of their choosing and they'll need to set the `wm.bar.status_config` to the absolute path of their config file in their XResources.  See the [upstream project documentation](https://docs.rs/i3status-rs/latest/i3status_rs/blocks/index.html) for further details.
