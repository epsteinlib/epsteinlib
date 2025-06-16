---
title: "1.5 Release Notes"
linkTitle: "1.5 Release Notes"
weight: 95
description: >
  Release notes for Regolith 1.5
---

Regolith R1.5 is a feature release which includes several improvements and optimizations. To summarize, Regolith 1.5 ships simpler workspace management, a Rofi-based Look switcher, and numerous internal optimizations and cleanup. Read below for more details.

## Features

| Feature                                 | Screenshot                                                                                                                                                                                                           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Next Free Workspace                     |                                                                                                                                                                                                                      | A typical part of managing workspaces in an i3-based desktop is moving to unused workspaces and then loading some applications. Before this feature, a user has to determine which unused workspace they prefer. This is done by scanning the list of existing used workspaces to determine an unused one. Now, the system can do this automatically. The {{< keys "super,`" >}} keybinding will move to the next free workspace. {{< keys "super,alt,`" >}} will move the focused window into the next free workspace.                                                                                          |
| View and Change Looks via Rofi          | {{< img "images/releases/regolith-select-look.png" "Select Rofi look" >}}                                                                                                                                            | Looks can be changed now via a Rofi dialog rather than having to configure the `Xresources` override via the command-line. To do this, use the keybinding {{< keys "super,alt,l" >}} and then select a look from the dialog.                                                                                                                                                                                                                                                                                                                                                                                     |
| GSettings Overrides                     |                                                                                                                                                                                                                      | Regolith now uses [GSettings overrides](https://help.gnome.org/admin/system-admin-guide/stable/overrides.html.en) to configure various GNOME settings for use with Regolith. In previous versions, settings were written globally to the user session from within the Regolith startup code. This could cause issues if the user works in multiple desktop environments. Now, the Regolith GNOME settings are defined in an override file that is only in effect while using a Regolith session. This allows switching between desktop environments without settings from Regolith impacting other environments. |
| New Looks                               | {{< img "images/releases/regolith-dracula.png" "Dracula Look" >}}<br>&nbsp;{{< img "images/releases/regolith-gruvbox.png" "Gruvbox Look" >}}<br>&nbsp;{{< img "images/releases/regolith-pop-os.png" "PopOs Look" >}} | Users have contributed some new Looks to Regolith: Dracula, Gruvbox, and PopOs. Each of these looks presents a distinctive color palate, typeface, and GTK theme.                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| i3-gaps upgraded to 4.18.2              |                                                                                                                                                                                                                      | See the [i3-gaps release notes](https://github.com/Airblader/i3/blob/a4a1a44275ea402b25d2d1365e1163e496024358/RELEASE-NOTES-4.18.2) for more information                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| More Refined Customizations             |                                                                                                                                                                                                                      | Numerous small changes allow more granular system customization, such as specifying the temperature unit, custom Compositor settings, and a more comprehensive way of changing i3 keybindings without having to copy the entire config file.                                                                                                                                                                                                                                                                                                                                                                     |
| More Desktop Environment Packages       |                                                                                                                                                                                                                      | The following packages can be installed in place of `regolith-desktop` for specific sets of packages based on user needs: `regolith-desktop-minimal`, `regolith-desktop-standard`, `regolith-desktop-mobile`, and `regolith-desktop-complete`                                                                                                                                                                                                                                                                                                                                                                    |
| New default compositor: Picom version 8 |                                                                                                                                                                                                                      | See the [Picom release notes](https://github.com/yshui/picom/releases) for more information on what has changed.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Remontoire upgraded to version 1.4      |                                                                                                                                                                                                                      | Includes better multi-monitor support and other bug fixes and enhancements.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Optional integration with **td-cli**    | {{< img "images/releases/regolith-td.png" "td-cli" >}}                                                                                                                                                               | Access a simple todo app via Rofi.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Documentation of development process    |                                                                                                                                                                                                                      | The [Regolith development process](https://regolith-linux.org/docs/policy-and-process/development/) is now better documented for greater transparency and to foster inclusion.                                                                                                                                                                                                                                                                                                                                                                                                                                             |

## Fixes

Have a look at [the R1.5 project page](https://github.com/orgs/regolith-linux/projects/12) for a list of bug fixes.

## Changelog Delta from Regolith 1.4.1 to Regolith 1.5

```
########################################
# Release Notes for dracula-gtk
########################################
dracula-gtk (1.0.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove unnecessary files


dracula-gtk (1.0-1) bionic; urgency=medium

  * [ Ken Gilmer ]
  * Packaging version add4f8c

########################################
# Release Notes for fonts-materialdesignicons-webfont
########################################
fonts-materialdesignicons-webfont (1.6.50-3regolith3) bionic; urgency=medium

  * Backporting to bionic for Regolith.


########################################
# Release Notes for gruvbox-gtk
########################################
gruvbox-gtk (1.0.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Rename root directory of theme to Gruvbox for consistency w/ other GTK themes.
  * Add gbp config for package management.


gruvbox-gtk (1.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * New upstream version 1.0


gruvbox-gtk (1.0-1) bionic; urgency=medium

  [ eximus ]
  * Initial commit
  * gruvbox theme


########################################
# Release Notes for i3-gaps-wm
########################################
i3-gaps-wm (4.18.2-1~regolith2) bionic; urgency=medium

  * Package source from upstream https://github.com/Airblader/i3/releases/tag/4.18.2


########################################
# Release Notes for i3ipc-python
########################################
i3ipc-python (2.1.1-1ubuntu1~ppa6) bionic; urgency=medium

  * Update build dependencies hoping to resolve packaing problems.



########################################
# Release Notes for i3xrocks
########################################
i3xrocks (1.3.4-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Version bump to match changelog.  Cleanup.


i3xrocks (1.3.3-1) bionic; urgency=medium

  [ Will Winder ]
  * Add optional default resource value.
  * Minor cleanup.
  * Free resource allocated by xcb_xrm_resource_get_string
  * Fix possible truncated resource value.

  [ Ken Gilmer ]
  * Add gbp config file


########################################
# Release Notes for picom
########################################
picom (8-1~1.gbp353272ubuntu1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove github files from debian branch.


picom (8-1~1.gbp353272) bionic; urgency=medium

  [ Ken Gilmer ]
  * New upstream version 8


########################################
# Release Notes for plano-theme
########################################
plano-theme (3.36-1-1regolith1) bionic; urgency=medium

  [ Ken Gilmer ]
  * New upstream version 3.36-1


########################################
# Release Notes for plymouth-theme-regolith
########################################
plymouth-theme-regolith (1.0.3-1) focal; urgency=medium

  * Tweaks to config files.


plymouth-theme-regolith (1.0.2-1) focal; urgency=medium

  * Ship grub file.


plymouth-theme-regolith (1.0.1-1) focal; urgency=medium

  * Add package hooks.



########################################
# Release Notes for pop-fonts
########################################
pop-fonts (1.0.3~1555617065~18.04~a86eb73) bionic; urgency=medium

  * Auto Build

########################################
# Release Notes for python3-i3ipc
########################################
python3-i3ipc (2.1.1-1ubuntu1~ppa7) bionic; urgency=medium

  * Changes to python version.


i3ipc-python (2.1.1-1ubuntu1~ppa6) bionic; urgency=medium

  * Update build dependencies hoping to resolve packaing problems.


i3ipc-python (2.1.1-1ubuntu1~ppa4) eoan; urgency=medium

  * Add python-xlib dependency.


i3ipc-python (2.1.1-1ubuntu1~ppa2) eoan; urgency=medium

  * Initial release from https://github.com/altdesktop/i3ipc-python/archive/v2.1.1.tar.gz.


########################################
# Release Notes for regolith-compositor-compton-glx
########################################
regolith-compositor-compton-glx (1.1.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Init script will kill pre-existing app instance before starting new. Part of fix for https://github.com/regolith-linux/regolith-desktop/issues/475.


regolith-compositor-compton-glx (1.0.10-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix typo in compton config file, found by @gservat.


regolith-compositor-compton-glx (1.0.9-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add xrender-sync-fence to handle issue https://github.com/regolith-linux/regolith-desktop/issues/116.


regolith-compositor-compton-glx (1.0.8-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Ship config file.


regolith-compositor-compton-glx (1.0.7-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.



########################################
# Release Notes for regolith-compositor-none
########################################
regolith-compositor-none (1.0.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.



########################################
# Release Notes for regolith-compositor-picom-glx
########################################
regolith-compositor-picom-glx (1.1.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix typo found in https://github.com/regolith-linux/regolith-compositor-picom-glx/issues/1.
  * Init script will kill pre-existing app instance before starting new. Part of fix for https://github.com/regolith-linux/regolith-desktop/issues/475.


regolith-compositor-picom-glx (1.0.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix typo found in
  https://github.com/regolith-linux/regolith-compositor-picom-glx/issues/1.


regolith-compositor-picom-glx (1.0.0-1) bionic; urgency=medium

  * Initial release


########################################
# Release Notes for regolith-compositor-xcompmgr
########################################
regolith-compositor-xcompmgr (1.2.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Init script will kill pre-existing app instance before starting new. Part of fix for https://github.com/regolith-linux/regolith-desktop/issues/475.


regolith-compositor-xcompmgr (1.1.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add ability to override xcompmgr defaults. Fixes https://github.com/regolith-linux/regolith-desktop/issues/382.


regolith-compositor-xcompmgr (1.0.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.



########################################
# Release Notes for regolith-default-settings
########################################
regolith-default-settings (1.0-1bionic1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add bionic specific gsettings overrides.


regolith-default-settings (1.0-1) focal; urgency=medium

  * Initial release, files moved from package regolith-gnome-flashback.


########################################
# Release Notes for regolith-desktop
########################################
regolith-desktop (2.78-1bionic) bionic; urgency=medium

  [ Ken Gilmer ]
  * Move from compton to picom as default compositor.


regolith-desktop (2.75-1bionic) bionic; urgency=medium

  [ Ken Gilmer ]
  * Change unclutter to regolith-unclutter-fixes to address https://github.com/regolith-linux/speed-ricer/issues/14.


regolith-desktop (2.72-1bionic1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Version bump for next DE release


regolith-desktop (2.71-1bionic1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Depend on regolith-default-settings to migrate to gsettings overrides.
  * Update branch name to conform w/ gbp.


regolith-desktop (2.68-1bionic1) bionic; urgency=medium

  * Package: Remove postinst to unblock upgrades (#432)


regolith-desktop (2.60-1bionic1) bionic; urgency=medium

  * Package: Remove unclutter-startup, if installed


regolith-desktop (2.59-1bionic1) bionic; urgency=medium

  * Bump version file for R1.4.1 release.



########################################
# Release Notes for regolith-ftue
########################################
regolith-ftue (1.0.11-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp config file.


########################################
# Release Notes for regolith-gnome-flashback
########################################
regolith-gnome-flashback (2.6.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove invalid component from desktop session file.
  * Cleanup


regolith-gnome-flashback (2.6.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Rename files to make session work.


regolith-gnome-flashback (2.6.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp config file.
  * Simply package, model directories in package and remove Makefile.
  * Update session files to refer to regolith over i3-gnome-flashback.  Experimental.


regolith-gnome-flashback (2.5.0-1) bionic; urgency=medium

  * Remove desktop files moved to regolith-default-settings.


regolith-gnome-flashback (2.4.19-1focal1) focal; urgency=medium

  * Fix wallpaper gsettings override.


regolith-gnome-flashback (2.4.18-1focal1) focal; urgency=medium

  * Remove non Regolith i3 config file load logic.


regolith-gnome-flashback (2.4.17-1focal1) focal; urgency=medium

  * Move to gsettings overrides for static dconf changes.
  * Change session name from i3-gnome to regolith.
  * Remove wm theme setting.



########################################
# Release Notes for regolith-i3-gaps-config
########################################
regolith-i3-gaps-config (2.8.0-1) bionic; urgency=medium

  [ Bence Ferdinandy ]
  * adding bindings for next_on_output and prev_on_output (#29)


regolith-i3-gaps-config (2.7.2-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add nm-applet to launch by default. Enables better network connectivity workflow.


regolith-i3-gaps-config (2.7.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Shorten look-selector to look for rofi menu.


regolith-i3-gaps-config (2.7.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add launcher for look selector.


regolith-i3-gaps-config (2.6.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Based on idea from @Ukaszg, extract all keybindings to Xres keys for easy overriding without i3 file copy.  Fixes https://github.com/regolith-linux/regolith-i3-gaps-config/issues/25.
  * Restart compositor upon i3 restart.  Part of fix for https://github.com/regolith-linux/regolith-desktop/issues/475.


regolith-i3-gaps-config (2.5.1-1) bionic; urgency=medium

  [ Moritz Heiber ]
  * Set tray output to primary by default (#24)

  [ Ken Gilmer ]
  * Allow the bar mode (dock|hide|invisible) to be specified by xres overrides.


regolith-i3-gaps-config (2.5.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add keybindings for 'next free workspace' feature.


regolith-i3-gaps-config (2.4.18-1) bionic; urgency=medium

  [ Will Winder ]
  * Add missing keybind for moving the workspace with hjkl. (#19)

  [ Jordan Cohen ]
  * shortcuts for explicit window orientation and hierarchical focusing (#18)


regolith-i3-gaps-config (2.4.17-1) bionic; urgency=medium

  [ Markus Bertheau ]
  * Allow customization of workspace_auto_back_and_forth (#20)


regolith-i3-gaps-config (2.4.16-1) bionic; urgency=medium

  [ ukaszg ]
  * configurable workspace keys (#17)


regolith-i3-gaps-config (2.4.15-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Cleanup from lintian warnings.


regolith-i3-gaps-config (2.4.14-1) bionic; urgency=medium

  * Take smart-gaps fix from @winder.
  * Comment cleanup.



########################################
# Release Notes for regolith-i3xrocks-config
########################################
regolith-i3xrocks-config (3.2.6-1ubuntu2) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix todo icon from updated material design web icons change


regolith-i3xrocks-config (3.2.5-1ubuntu1) bionic; urgency=medium

  [ Rafael Lacerda ]
  * Fix crash when focused window name is None (#81)

  [ Peter Mahnke ]
  * Remove decimal places from battery number (#80)


regolith-i3xrocks-config (3.2.4-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Temp status indicator can be configured with sensor name via Xres value i3xrocks.temp.key.


regolith-i3xrocks-config (3.2.3-1) bionic; urgency=medium

  [ Pedro Kostelec ]
  * [IMPR] Show aggregated battery percentage (#78)

  [ Ken Gilmer ]
  * Fix issue found in livecd enviroment w memory script

  [ Pedro Kostelec ]
  * [FIX] Fixes dissapearing indicator when value is 0 (#79)


regolith-i3xrocks-config (3.2.2-1) bionic; urgency=medium

  [ Peter Mahnke ]
  * Make curl work for special characters (#76)

  [ Ken Gilmer ]
  * Remove debugging line mistakenly added by me in earlier fix.


regolith-i3xrocks-config (3.2.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix cpu-usage script


regolith-i3xrocks-config (3.2.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add ability to specify temperature unit via xres, and fallback to country locale if unspecified. Based on https://github.com/regolith-linux/regolith-desktop/issues/372. (#72)

  [ Joel Hoisko ]
  * Add option to define a specific media player (#63)


regolith-i3xrocks-config (3.1.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add ability to override click commands via Xresources. Fixes https://github.com/regolith-linux/regolith-desktop/issues/352.


regolith-i3xrocks-config (3.1.2-1) bionic; urgency=medium

  [ Leo Palmer Sunmo ]
  * Polish up media-player text splitting and pango formatting (#71)


regolith-i3xrocks-config (3.1.1-1) bionic; urgency=medium

  [ Vipul Mathur ]
  * Use `MemAvailable` to indicate memory usage (#70)


regolith-i3xrocks-config (3.1.0-1) bionic; urgency=medium

  [ Bence Ferdinandy ]
  * fixed an error with name somehow being None
  * adding next-workspace


regolith-i3xrocks-config (3.0.34-1) bionic; urgency=medium

  [ Marten Lienen ]
  * Exit gracefully if the focused window does not have a name (#68)


regolith-i3xrocks-config (3.0.33-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * regress debhelper-compat version due to build failure in launchpad.


regolith-i3xrocks-config (3.0.32-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Downgrade debhelper version to appease launchpad.
  * Remove compat file in attempt to appease debuild.


regolith-i3xrocks-config (3.0.31-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Downgrade debhelper version to appease launchpad.


regolith-i3xrocks-config (3.0.30-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.


regolith-i3xrocks-config (3.0.29-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Package: @cheginit's new todo status indicator based on td-cli.


regolith-i3xrocks-config (3.0.28-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Package: @buck10's major rework of keyboard-layout.


regolith-i3xrocks-config (3.0.27-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Package: Removed verbose output from net-traffic script (#62)


regolith-i3xrocks-config (3.0.26-2) bionic; urgency=medium

  [ Ken Gilmer ]
  * Update i3xrocks-focused-window-name to use consistent python library name.


regolith-i3xrocks-config (3.0.26-1) bionic; urgency=medium

  [ Taher Chegini ]
  * Add todo blocklet (#54)

  [ Ken Gilmer ]
  * Initial upstream branch.
  * New upstream version 3.0.25


regolith-i3xrocks-config (3.0.25-2) bionic; urgency=medium

  * Trying w/ debhelper v11 based on build error.


regolith-i3xrocks-config (3.0.25-1) bionic; urgency=medium

  * Package changes by @moritzheiber:
  ** Fix button clicks blocklets
  ** Also fix network bandwidth calculation/shellcheck errors in net-traffic script


regolith-i3xrocks-config (3.0.24-1) bionic; urgency=medium

  * Package @moritzheiber's fixes for disk-capacity and other blocks.
  * Package @cjohnston1158 additions to net-traffic.


regolith-i3xrocks-config (3.0.23-1) bionic; urgency=medium

  * Add info indicator missed when moving to conf.d config.


regolith-i3xrocks-config (3.0.22-2) bionic; urgency=medium

  * Fixes to weather status for ambiguous city names and cities of multiple
    words.


########################################
# Release Notes for regolith-rofi-config
########################################
regolith-rofi-config (1.3.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix look-selector


regolith-rofi-config (1.3.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add script to view and select Looks.  Minor cleanup.  Addresses https://github.com/orgs/regolith-linux/projects/12.


regolith-rofi-config (1.2.6-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp config file.



########################################
# Release Notes for regolith-rofication
########################################
regolith-rofication (1.2.3-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Remove scripts that do not apply to user sessions


regolith-rofication (1.2.2-2) focal; urgency=medium

  * Add git-buildpackage config


########################################
# Release Notes for regolith-styles
########################################
regolith-styles (2.6.13-1ubuntu2) bionic; urgency=medium

  [ Ken Gilmer ]
  * More tweaks to icons for material design based looks


regolith-styles (2.6.12-1ubuntu1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Update material design icon-based looks to use fonts-materialdesignicons-webfont 1.6 mappings


regolith-styles (2.6.11-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add regolith.look Xres entry for currently set look.


regolith-styles (2.6.10-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * If Regolith user root Xres file exists, move old copy to tmp instead of aborting.  Sort looks.  Misc cleanup.


regolith-styles (2.6.9-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Check a flag before executing gnome-terminal profile script.  Fixes https://github.com/regolith-linux/regolith-desktop/issues/471.


regolith-styles (2.6.8-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Update term colors for Dracula for more variance and remove conflict for color 11.  Addresses https://github.com/regolith-linux/regolith-styles/issues/32.


regolith-styles (2.6.7-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Change selected item background color to addresss https://github.com/regolith-linux/regolith-styles/issues/23.


regolith-styles (2.6.6-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add next-workspace glyph for missing looks.  Change sound, mute, ram glyphs for fontawesome looks.


regolith-styles (2.6.5-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix font dependency for pop-os.  Cleanup control file.


regolith-styles (2.6.4-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add package for regolith-look-pop-os.


regolith-styles (2.6.3-1) bionic; urgency=medium

  [ Aditya Wagh ]
  * Create Pop!_OS Style (#29)

  [ Ken Gilmer ]
  * Update gtk theme name for consistency.
  * Add packaging data for pop-os look.


regolith-styles (2.6.2-1) bionic; urgency=medium

  * Specify glyph for next-workspace status indicator for fontawesome and materialdesign.


regolith-styles (2.6.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Package the dracula look.


regolith-styles (2.6.0-1) bionic; urgency=medium

  [ Drew DeVore ]
  * Create Dracula Theme (#30)

  [ Julius Francisco ]
  * Remove user specific path (#31)


regolith-styles (2.4.27-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix typo in package file.


regolith-styles (2.4.26-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix gruvbox missing files on install.


regolith-styles (2.4.25-1) bionic; urgency=medium

  [ Ekans ]
  * Create gruvbox style (#25)

  [ Ken Gilmer ]
  * Add packaging for gruvbox.


regolith-styles (2.4.24-1) bionic; urgency=medium

  [ Chua Hou ]
  * Balance workspace name by adding space at end (#20)


regolith-styles (2.4.23-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.



########################################
# Release Notes for regolith-system
########################################
regolith-system (1.4.0.4-1) focal; urgency=medium

  [ Ken Gilmer ]
  * Add gbp conf file, cleanup.



########################################
# Release Notes for regolith-unclutter-xfixes
########################################
regolith-unclutter-xfixes (1.5-2) bionic; urgency=medium

  [ Ken Gilmer ]
  * Fix linter warning
  * Revert "Remove master files"


regolith-unclutter-xfixes (1.5-1) disco; urgency=medium

  [ Ken Gilmer ]
  * Initial release


########################################
# Release Notes for remontoire
########################################
remontoire (1.4.0-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Specify font for windows glyph in stylesheet.
  * Increase buffer to 128Kb to handle larger i3 config files
  * Add feature to consume for stdin and specify arbitrary string prefix for more config file format support.
  * Fix window positioning in multi-monitor setups.


########################################
# Release Notes for td-cli
########################################
td-cli (1.2.4-7) bionic; urgency=medium

  * Initial release


########################################
# Release Notes for xrescat
########################################
xrescat (1.2.1-1) bionic; urgency=medium

  [ Ken Gilmer ]
  * Cleanup lintian warnings.
  * Version bump to appease launchpad


xrescat (1.2.0-1) bionic; urgency=medium

  [ Jeremy Coghill ]
  * added logic to exclude surrounding double quotes from returned strings to allow for the use of leading and trailing spaces in resources. updated version number. updated man page

  [ Ken Gilmer ]
  * Update README for new double quote handling.  Add gbp config. Cleanup.

```
