---
title: "Manage a To-Do List"
linkTitle: "Manage a To-Do List"
description: >
  Create and manage a to-do list from status bar
---

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

[`td-cli`](https://github.com/darrikonn/td-cli) is a command line todo manager written in Python3 for organizing todo lists across multiple projects.

Regolith includes a status bar blocklet for accessing `td-cli` directly from the bar. It can be installed as follows:

```bash
sudo apt install i3xrocks-todo
```

After installing `td-cli` you can access it via command line using `td`. The full documentation of `td-cli`'s API can be found [online at their repository](https://github.com/darrikonn/td-cli/blob/master/API.md).
For example, let's create a todo list with a group name of Regolith, set it as the default project, and then add a new task to this group:

```bash
td ag regolith
td g regolith p
td a "Write a Howto for td-cli"
```

Now, if you click on the `td-cli` blocklet (on the left of the time blocklet in the figure below) a floating terminal opens up with an interactive instance of `td-cli`:

{{< img "images/tdcli.png" "The td-cli indicator on the i3 bar">}}

{{< img "images/tdcli_window.png" "The td-cli pop-out window">}}

You can set this blocklet to only show number of uncompleted tasks in the bar by setting the following key in your `Xresource` file, for example, `~/.config/regolith3/Xresources`:

```yaml {filename="~/.config/regolith3/Xresources"}
i3xrocks.todo.uncompleted: true
```

By default, this key is set to `false`. Moreover, you can add a keybinding for accessing `td-cli` via a floating terminal by adding the following line to your `i3` config file, for example, `~/.config/regolith3/i3/config`:

```text {filename="~/.config/regolith3/i3/config"}
bindsym $mod+Ctrl+d exec --no-startup-id "/usr/bin/gnome-terminal --class=floating_window -e 'td --interactive'"
```

Alternatively, you add these lines to your `i3` config file, first:

```text {filename="~/.config/regolith3/i3/config"}
set_from_resource $wm.bindsym.1 wm.bindsym.1 :
set_from_resource $wm.bindsym.program.1 wm.bindsym.program.1 :
bindsym $mod+$wm.bindsym.1 exec --no-startup-id "$wm.bindsym.program.1"
```

Then, add the following lines to your `Xresources` file:

```yaml {filename="~/.config/regolith3/Xresources"}
wm.bindsym.1: Ctrl+d
wm.bindsym.program.1: /usr/bin/gnome-terminal --class=floating_window -e 'td --interactive'
```
