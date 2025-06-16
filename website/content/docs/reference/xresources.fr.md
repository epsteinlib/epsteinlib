---
title: "Xresources"
linkTitle: "Xresources"
weight: 7
description: >
  Toutes les clés `Xresources` disponibles pour Regolith
---

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

## Main `Xresources` keys

| `Xresources` Key                                   | Description                                                                                                                   |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| gtk.theme_name                                     | GTK Theme ("Legacy Applications" in `gnome-tweaks`) |
| gtk.icon_theme_name                                | GTK Icon Theme                                                                                                                |
| gnome.terminal.audible-bell                        | Enables audible bell in gnome-terminal (`true / false`)                                                                       |
| gtk.document_font_name                             | GTK Document Font Name                                                                                                        |
| gtk.font_name                                      | General purpose font for GTK apps and the desktop                                                                             |
| gtk.monospace_font_name                            | Font for monospace areas such as i3bar and the terminal                                                                       |
| gnome.terminal.background-transparency-percent     | Background transparency value for gnome-terminal (`0 - 100`). See also `gnome.terminal.use-transparent-background`.           |
| gnome.terminal.font                                | Font definition for gnome-terminal                                                                                            |
| gnome.terminal.scrollbar                           | Flag to enable/disable scrollbar (`never / always`)                                                                           |
| gnome.terminal.use-transparent-background          | Enables transparent background in gnome-terminal (`true / false`). See also `gnome.terminal.background-transparency-percent`. |
| gnome.wallpaper                                    | GNOME Wallpaper (can be overridden in Settings)                                                                               |
| gnome.wm.theme                                     | GNOME Window Manager Theme (unused)                                                                                           |
| wm.alt                                          | Key mapping for Alt key. See the [i3 User Guide](https://i3wm.org/docs/userguide.html#_using_i3) for details.                 |
| wm.bar.background.color                         | Bar Background Color                                                                                                          |
| wm.bar.position                                 | Bar Screen Position                                                                                                           |
| wm.bar.separator.color                          | Color of separator between blocks                                                                                             |
| wm.bar.statusline.color                         |                                                                                                                               |
| wm.bar.trayoutput                               | Enable or disable the [app tray](https://i3wm.org/docs/userguide.html#_tray_output) in i3bar                                  |
| wm.bar.workspace.active.background.color        |                                                                                                                               |
| wm.bar.workspace.active.border.color            |                                                                                                                               |
| wm.bar.workspace.active.text.color              |                                                                                                                               |
| wm.bar.workspace.focused.background.color       |                                                                                                                               |
| wm.bar.workspace.focused.border.color           |                                                                                                                               |
| wm.bar.workspace.focused.text.color             |                                                                                                                               |
| wm.bar.workspace.inactive.background.color      |                                                                                                                               |
| wm.bar.workspace.inactive.border.color          |                                                                                                                               |
| wm.bar.workspace.inactive.text.color            |                                                                                                                               |
| wm.bar.workspace.urgent.background.color        |                                                                                                                               |
| wm.bar.workspace.urgent.border.color            |                                                                                                                               |
| wm.bar.workspace.urgent.text.color              |                                                                                                                               |
| wm.client.focused.color.background              |                                                                                                                               |
| wm.client.focused.color.border                  |                                                                                                                               |
| wm.client.focused.color.child_border            |                                                                                                                               |
| wm.client.focused.color.indicator               |                                                                                                                               |
| wm.client.focused.color.text                    |                                                                                                                               |
| wm.client.focused_inactive.color.background     |                                                                                                                               |
| wm.client.focused_inactive.color.border         |                                                                                                                               |
| wm.client.focused_inactive.color.child_border   |                                                                                                                               |
| wm.client.focused_inactive.color.indicator      |                                                                                                                               |
| wm.client.focused_inactive.color.text           |                                                                                                                               |
| wm.client.unfocused.color.background            |                                                                                                                               |
| wm.client.unfocused.color.border                |                                                                                                                               |
| wm.client.unfocused.color.child_border          |                                                                                                                               |
| wm.client.unfocused.color.indicator             |                                                                                                                               |
| wm.client.unfocused.color.text                  |                                                                                                                               |
| wm.client.urgent.color.background               |                                                                                                                               |
| wm.client.urgent.color.border                   |                                                                                                                               |
| wm.client.urgent.color.child_border             |                                                                                                                               |
| wm.client.urgent.color.indicator                |                                                                                                                               |
| wm.client.urgent.color.text                     |                                                                                                                               |
| wm.floatingwindow.border.size                   | Window border size in pixels for floating windows                                                                             |
| wm.font                                         | Window manager font                                                                                                           |
| wm.gaps.inner.size                              | Default gap size between windows in pixels                                                                                    |
| wm.gaps.outer.size                              | Default gap size at screen edge in pixels                                                                                     |
| wm.mod                                          | Key mapping for Super key. See the [i3 User Guide](https://i3wm.org/docs/userguide.html#_using_i3) for details.               |
| wm.window.border.size                           | Default border size in pixels                                                                                                 |
| wm.workspace.01.name .. wm.workspace.19.name | Workspace Labels                                                                                                              |
| wm.workspace.01.key .. wm.workspace.10.key   | Workspace Shortucts                                                                                                           |
| i3xrocks.action.volume.left                        | Command to execute when volume indicator is left clicked (default: none)                                                      |
| i3xrocks.action.volume.middle                      | Command to execute when volume indicator is middle clicked (default: none)                                                    |
| i3xrocks.action.volume.right                       | Command to execute when volume indicator is right clicked (default: mute through amixer)                                      |
| i3xrocks.action.volume.scrollup                    | Command to execute when scroll wheel up is used on volume indicator (default: increase volume through amixer)                 |
| i3xrocks.action.volume.scrolldn                    | Command to execute when scroll wheel down is used on volume indicator (default: decrease volume through amixer)               |
| i3xrocks.error.color                               | Color for bar module in error mode                                                                                            |
| i3xrocks.label.battery0                            | Bar icon glyph                                                                                                                |
| i3xrocks.label.battery100                          | Bar icon glyph                                                                                                                |
| i3xrocks.label.battery20                           | Bar icon glyph                                                                                                                |
| i3xrocks.label.battery50                           | Bar icon glyph                                                                                                                |
| i3xrocks.label.battery80                           | Bar icon glyph                                                                                                                |
| i3xrocks.label.batteryac                           | Bar icon glyph                                                                                                                |
| i3xrocks.label.batterycharging                     | Bar icon glyph                                                                                                                |
| i3xrocks.label.batteryunknown                      | Bar icon glyph                                                                                                                |
| i3xrocks.label.capslock                            | Bar icon glyph                                                                                                                |
| i3xrocks.label.color                               | Bar icon glyph                                                                                                                |
| i3xrocks.label.cpu                                 | Bar icon glyph                                                                                                                |
| i3xrocks.label.disk                                | Bar icon glyph                                                                                                                |
| i3xrocks.label.dn                                  | Bar icon glyph                                                                                                                |
| i3xrocks.label.help                                | Bar icon glyph                                                                                                                |
| i3xrocks.label.mediapaused                         | Bar icon glyph                                                                                                                |
| i3xrocks.label.mediaplaying                        | Bar icon glyph                                                                                                                |
| i3xrocks.label.memory                              | Bar icon glyph                                                                                                                |
| i3xrocks.label.mute                                | Bar icon glyph                                                                                                                |
| i3xrocks.label.notify.error                        | Bar icon glyph                                                                                                                |
| i3xrocks.label.notify.none                         | Bar icon glyph                                                                                                                |
| i3xrocks.label.notify.some                         | Bar icon glyph                                                                                                                |
| i3xrocks.label.numlock                             | Bar icon glyph                                                                                                                |
| i3xrocks.label.off                                 | Bar icon glyph                                                                                                                |
| i3xrocks.label.on                                  | Bar icon glyph                                                                                                                |
| i3xrocks.label.sound                               | Bar icon glyph                                                                                                                |
| i3xrocks.label.thermometer                         | Bar icon glyph                                                                                                                |
| i3xrocks.label.time                                | Bar icon glyph                                                                                                                |
| i3xrocks.label.up                                  | Bar icon glyph                                                                                                                |
| i3xrocks.label.wifi                                | Bar icon glyph                                                                                                                |
| i3xrocks.label.wired                               | Bar icon glyph                                                                                                                |
| i3xrocks.nominal                                   | Default bar icon color                                                                                                        |
| i3xrocks.value.color                               | Color of numerical value in bar module                                                                                        |
| i3xrocks.value.font                                | Font for numerical value in bar module                                                                                        |
| i3xrocks.warning                                   | Warning color for bar module                                                                                                  |

## Optional `Xresources` keys

The following `Xresources` keys are undefined by default but can be used by users to change behavior as desired.

| `Xresources` Key                    | Description                                                                                                                                                                                         |
| ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| wm.gaps.focus_follows_mouse      | The window under the mouse has focus if enabled                                                                                                                                                     |
| wm.program.launcher.app          | Command i3 executes from the app launcher keybinding                                                                                                                                                |
| wm.program.launcher.cmd          | Command i3 executes from the command launcher keybinding                                                                                                                                            |
| wm.program.launcher.window       | Command i3 executes for the window navigation keybinding                                                                                                                                            |
| wm.program.help                  | Command to launch from the help keybinding                                                                                                                                                          |
| wm.program.file_search           | Command to launch from the file search keybinding                                                                                                                                                   |
| wm.program.refresh_ui            | Command to launch from the session refresh keybinding                                                                                                                                               |
| wm.program.logout                | Command to launch from the logout keybinding                                                                                                                                                        |
| wm.program.reboot                | Command to launch from the reboot keybinding                                                                                                                                                        |
| wm.program.shutdown              | Command to launch from the shutdown keybinding                                                                                                                                                      |
| wm.program.lock                  | Command to launch from the lock screen keybinding                                                                                                                                                   |
| wm.program.sleep                 | Command to launch from the sleep keybinding                                                                                                                                                         |
| wm.program.settings              | Command to launch from the Control Panel keybinding                                                                                                                                                 |
| wm.program.display               | Command to launch from the Display Settings keybinding                                                                                                                                              |
| wm.program.wifi                  | Command to launch from the Wifi Settings keybinding                                                                                                                                                 |
| wm.program.bluetooth             | Command to launch from the Bluetooth Settings keybinding                                                                                                                                            |
| wm.program.files                 | Command to launch from the file browser keybinding                                                                                                                                                  |
| wm.program.notification_ui       | Command to launch from the notification UI keybinding                                                                                                                                               |
| wm.program.1...wm.program.3   | Optional programs to start with i3                                                                                                                                                                  |
| i3xrocks.media-player.player-args   | Optional arguments used by the `media-player` i3xrocks block. These arguments are fed straight into `playerctl`, eg. setting this value to `-p spotify` only makes the bar display/control spotify. |
| i3xrocks.action.mic.left            | Optional command to run on left button click of microphone status indicator.                                                                                                                        |
| i3xrocks.action.mic.middle          | Optional command to run on middle button click of microphone status indicator.                                                                                                                      |
| i3xrocks.action.mic.right           | Override command to run on right button click of microphone status indicator.                                                                                                                       |
| i3xrocks.action.mic.scrollup        | Override command to run on scroll up of microphone status indicator.                                                                                                                                |
| i3xrocks.action.mic.scrolldn        | Override command to run on scroll down of microphone status indicator.                                                                                                                              |
| i3xrocks.action.volume.left         | Optional command to run on left button click of microphone status indicator.                                                                                                                        |
| i3xrocks.action.volume.middle       | Optional command to run on middle button click of microphone status indicator.                                                                                                                      |
| i3xrocks.action.volume.right        | Override command to run on right button click of microphone status indicator.                                                                                                                       |
| i3xrocks.action.volume.scrollup     | Override command to run on scroll up of microphone status indicator.                                                                                                                                |
| i3xrocks.action.volume.scrolldn     | Override command to run on scroll down of microphone status indicator.                                                                                                                              |
| i3xrocks.action.disk-capacity       | Override command to run on disk capacity status indicator.                                                                                                                                          |
| i3xrocks.action.cpu-usage           | Override command to run on cpu usage status indicator.                                                                                                                                              |
| i3xrocks.action.memory              | Override command to run on memory status indicator.                                                                                                                                                 |
| i3xrocks.action.wifi                | Override command to run on wifi status indicator.                                                                                                                                                   |
| i3xrocks.action.weather             | Override command to run on weather status indicator.                                                                                                                                                |
| i3xrocks.action.media-player.left   | Optional command to run on left button click of media player status indicator.                                                                                                                      |
| i3xrocks.action.media-player.middle | Optional command to run on middle button click of media player status indicator.                                                                                                                    |
| i3xrocks.action.media-player.right  | Override command to run on right button click of media player status indicator.                                                                                                                     |
| i3xrocks.action.time                | Override command to run on date time status indicator.                                                                                                                                              |
| i3xrocks.action.battery             | Override command to run on battery status indicator.                                                                                                                                                |
