---
title: "Basic Usage"
linkTitle: "Basic Usage"
weight: 3
description: >
  Learn the most important parts of Regolith first.
---

## Launch an Application

After the first login, the screen will contain a bar at the bottom and some keyboard shortcuts in a dialog.  This dialog, [`ilia`](https://github.com/regolith-linux/ilia), provides assistance for finding keybindings for various actions in the desktop environment.  It may be dismissed with the `esc` key and invoked again with {{< keys "super,shift,?" >}}. To launch an application, press {{< keys "super,space" >}} and a similar dialog will load with the applications installed on the system. If you begin typing the name of the application you wish to launch, the list will filter down to only those that match. You can also scroll through the list with the arrow keys.

![](/images/v-tour/regolith-ilia-apps.png)

## Navigate to a Window

To focus on a running application, the {{< keys "super,ctrl,space" >}} keybinding can be used to show a dialog of all running windows. Simply select the desired window and the system will move to the window's containing workspace and set its focus.

![](/images/v-tour/regolith-ilia-windows.png)

## Create a new Workspace

Workspaces (also known as "virtual desktops" in some desktop environments) are a primary way by which windows can be managed. Each workspace provides a screen for an independent set of windows. To create or navigate to an existing workspace, specify it's index with keys {{< keys "super,0-9" >}}.  

## Find and Launch a File

Searching your local computer for a specific file can be achieved in a variety of ways. Regolith provides a file search dialog via {{< keys "super,alt,space" >}} to quickly seek and optionally load a file with its default application. Pressing enter on a selection will cause the file to be invoked with the `xdg-open` command.

![](/images/v-tour/regolith-ilia-files.png)

## Terminal

The terminal is a special application that has its own keybinding, {{< keys "super,enter" >}}.

## Browser

The preferred browser is also special and can be launched with {{< keys "super,shift,enter" >}}.

## Quit an Application

{{< keys "super,shift,q" >}} will quit the focused application. For misbehaving apps, {{< keys "super,alt,q" >}} will force quit (kill) the focused application, but note that any unsaved data may be lost.

## Switch Workspaces

A primary way of managing how applications occupy the screen is by the usage of workspaces. A workspace is in essence a group of zero or more applications that all present windows on the same virtual display. Users can switch workspaces to gain more space for applications to occupy. Windows can be moved across workspaces. Workspaces can be tied to physical display or be moved to other displays as needed. To switch workspaces, hold down the {{< keys "super" >}} key and press one of the number keys. This will navigate to that workspace. If any windows exist in a given workspace, it will show in the bar as being active. Thus, by referring to the bar the user can know in which workspaces active applications live.

Workspaces can also be paged forward and backward via {{< keys "super,tab" >}} and {{< keys "super,shift,tab" >}}.

## Configuring the System

The Settings app (also known as `regolith-control-center`) is the primary interface for common configuration tasks, such as configuring WiFi and Bluetooth, configuring app associations, and setting up users. The app can be launched via the {{< keys "super,c" >}} keybinding or from the launcher via the `Regolith Settings` app.

![](/images/v-tour/regolith-gnome-settings.png)

## View Keybindings

All keybindings used by the Regolith desktop can be viewed, sorted, and executed by [`ilia`](https://github.com/regolith-linux/ilia).  This dialog is launched with {{< keys "super,shift,?" >}}.

![](/images/v-tour/regolith-ilia-keybindings.png)

## Manage Notifications

Regolith's desktop notification system is designed to prevent pop-ups that may distract the user from their work.  Notifications can be viewed and managed in [`ilia`](https://github.com/regolith-linux/ilia) via {{< keys "super,n" >}}.  Specific keybindings such as deleting one or groups of notifications may be found on the `help` tab of the dialog.  For users wishing for more visibility of current notifications, the status indicator `i3xrocks-rofication` may be installed via the system package manager. This will show the number of critical and normal unread notifications on the bar. 
