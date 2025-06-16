---
title: "First Launch"
linkTitle: "First Launch"
weight: 2
description: >
  Important details for the first time you start Regolith
---

Regolith's default interface lacks any visual hint as to how to launch applications or perform other common tasks. This can make the first time experience a bit confusing. However using the Regolith interface is relatively simple by only committing the following two keybindings to memory:

| Action | Keybinding | Purpose |
|--------|------------|---------|
| Keybinding Viewer | {{< keys "super,shift,?" >}} | View and execute available keybindings |
| App Launcher | {{< keys "super,space" >}} | Launch a desktop application |

Entering one of these keybindings will cause a launcher dialog to pop up on the focused monitor. Both of these dialogs have the same general interface.

The interface, a program called [`ilia`](https://github.com/regolith-linux/ilia), has a text entry widget at the top, followed by a list selection widget.  The general usage pattern for these dialogs is to select an item from the displayed list.  Begin typing to filter the selections, and use the arrow keys to make a selection.  Pressing {{< keys "enter" >}} will cause the selected item to be activated.  For keybindings, this means a keybinding will be manually executed.  For apps, the selected application will be launched.

## The Super Key

Most keybindings you'll be using with Regolith start with the {{< keys "super" >}} key. This key typically says "WIN" or has a MS Windows logo on it. On a Mac it is the Command {{< keys "⌘" >}} key. Going forward, this key will simply be referred to as {{< keys "super" >}} but keep in mind that you can change which physical button {{< keys "super" >}} maps to in your configuration if you choose to.

## Keybinding Viewer

![](/images/v-tour/regolith-ilia-keybinding-window.png)

The keybinding viewer's purpose is to make finding and learning the Regolith keybindings quick and easy.  A reference to the keybindings is available on this site, however it's faster to simply bring up the dialog and search for the action you are looking for.  Once found,
the keybinding can be executed either with the keyboard or by selecting the item in the dialog.

## Application Launcher

![](/images/v-tour/regolith-ilia-apps-window.png)

A fundamental task of any desktop environment is to enable users to launch applications.  In Regolith this is done via the app launcher dialog.  As with keybindings, the dialog first shows a list of all available desktop applications.  Begin typing to filter the selection and use the arrow keys to specify a specific selection.  Press {{< keys "enter" >}} to launch the selected application.

## Next Steps

Have a look at the [basics](../basics) page to learn more about the Regolith interface.
