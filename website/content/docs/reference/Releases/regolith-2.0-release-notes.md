---
title: "2.0 Release Notes"
linkTitle: "2.0 Release Notes"
weight: 90
description: >
  Release notes for Regolith 2.0
---

Regolith 2.0 is a major release with some breaking changes to looks and i3 configuration, and includes several improvements and optimizations.

## New Features

* Adding Debian Bullseye to the list of supported OS's including Ubuntu Focal and Jammy.
* The desktop UI component `rofi` was replaced with desktop executor designed specifically for Regolith: [`ilia`](https://github.com/regolith-linux/ilia).
* The desktop UI component `gnome-control-center` was replaced with a settings app designed specifically for Regolith: [`regolith-control-center`](https://github.com/regolith-linux/regolith-control-center).
* The Regolith i3 config file has been destructured into [partials](https://github.com/regolith-linux/regolith-i3-config/tree/master/partials) and packaged to provide better configuration management.
* The Regolith [Look XResource scheme](https://github.com/regolith-linux/regolith-look-default/blob/main/usr/share/regolith-look/default/root) has been updated to be simpler, more flexible, and consistent. 
* The keybindings viewer `remontoire` was replaced by `ilia` to provide better support for keyboard-driven workflow and easier searching.
* Regolith supports both `i3` and `i3-gaps` based on user preference.
* Fonts from the open source project [Nerd Fonts](https://nerdfonts.com) have been packaged for Regolith, simplifying font management.
* A command-line utility for managing notifications, `lago` is available.
* A command-line utility for collecting Regolith diagnostic info, `regolith-diagnostic` is available.

## Usability Fixes

* `Alt` and `Super` can now be swapped easily with an Xresources override.
* Alternate languages can be selected from in the desktop executor.
* `regolith-control-center` prevents some pop-ups from being hidden and locking UI

## Known Issues

* No localization for `ilia` 

## Changelog Delta from Regolith 2.0 Beta 3 to General Release

### Changes in `regolith-session`:
```
1061f67 Re-apply user Xres after look is applied. Resolves https://github.com/regolith-linux/regolith-desktop/issues/644
```

### Changes in `regolith-i3-config`:
```
e25d88c Add missing newlines at end of config files from @buck10
da5615f Add conflicts relationship between ilia and rofi i3 config
244547e Add descriptions for bin packages
a11ff93 Allow look-selector icon to be overridden. Specify better default icon
```

### Changes in `regolith-look-default`:
```
fa7fb24 Update ilia CSS style
659d545 Add docs to Xresources
```

### Changes in `regolith-look-extra`:
```
d7b784b Update ilia CSS style
bc88936 Fix look descriptions
fd84e99 Fix incorrect path to ilia css
```

### Changes in `regolith-i3xrocks-config`:
```
b35690d fixed xml escaping (#127)
```

## Changes in `ilia`:
```
fa548a8 Remove hard code style. Squelch right-click context menu to avoid hang.
9ed0020 Extract hard coded style to CSS
16d088b Add rofi reference
05dcaee Fix code sections
018deb4 Updates to readme for release
```

## Changes in `regolith-control-center`:
```
e85c3ced1 Fix crash on power panel (jammy only)
6a584746b Update regolith logo from final revision
```


## Changelog Delta from Regolith 2.0 Beta 2 to Beta 3

### Changes in `regolith-session`:
```
87919c5 Revert strategy that creates merged Xresource. Simplify by just loading one override or using default look. Looks become top-down definitions of UI, and bottom up extension model deprecated.
975c542 Cleanup
```

### Changes in `regolith-i3-config`:
```
4d0d5b2 polish look selector ui
```

### Changes in `regolith-look-default`:
```
557d76d Add lockscreen wallpaper functionality
```

### Changes in `regolith-look-extra`:
```
8c211ba Package changes
b8188d3 Visual tweaks to nevil
```

### Changes in `regolith-ftue`:
```
ce16a93 Use stylesheet for proper font rendering on ftue help
```

### Changes in `ilia`:
```
4512ba1 Add separator between category and action in keybindings view
b8373fd Exit if super pressed to prevent blocking i3
85f7041 Fix warning in socket close
a870a2b Better icon for window page
647b1bb Add missing dependency
e2c19c1 Update textlist to render icons
fb59b13 Plumb cmdline args to each page via DialogPage. Add title and icon parameters for TextList customization
b83fbcb Update dialog resize keys to use alt to allow keyboard entry. Specify max bounds of window size when possible.
```

### Changes in `regolith-control-center`:
```
3eb9dea29 Re-add some install files for shortcuts and gschema
```

## Changelog Delta from Regolith 2.0 Beta 1 to Beta 2

* New light look: regolith-look-nevil
* Updated wallpaper for `regolith-look-nord` and `regolith-look-gruvbox` by @GMFranceschini
* Ilia now has multi-tab mode
* Ilia launches app from notification

### Changes in `regolith-session`:
```
f57f66b Add script to collect diagnostic info for troubleshooting
4494e68 Fix broken refresh logic after Xres file change introduced in previous version
cdaa330 Fix regolith look update logic to write to Regolith-specific Xresource file
```

### Changes in `regolith-i3-config`:
```
406b428 Add docs about config partials
521016c Add OR dependency on i3-gaps-wm to resolve packaging issue found on Bullseye
a6b6fea Load all pages in ilia when users load the keybindings dialog
0f6a018 Add networkmanager, next-workspace, and unclutter as soft dependencies
7921136 Require i3-gaps version with config partials in regolith 2
a2390bd Require i3 version with config partials in regolith 2
2711991 Fix typo of new mapping for Xres mono font to i3 bar font
1238b4a Fix typo of new mapping for Xres font to i3 font
```

### Changes in `whitesur-gtk-theme`:
```
cc591d7 Initial commit. Built from release 2022-02-21
```

### Changes in `regolith-look-default`:
```
fcb288e Default terminal profile update to true
c849879 make the check for invalid path to wallpaper more robust
5c2c402 removed left-over bracket
ca61239 fixed regression in #2 + fixed other possible bug
ad6691b Merge pull request #2 from georglauterbach/main
c9df3c7 updated `default_loader.sh` | style, formatting & refactoring
```

### Changes in `regolith-look-extra`:
```
2e700f0 Use new wallpapers from GMFranceschini in nord and gruvbox looks
a2ca36b Fix dependency
5c7f3ed Tweaks and fixes to nevil
ec1a158 Add nevil look
429863c Add updated solarized-dark look
3511556 Give rounded corners to Nord w/ custom picom config.  Misc tweaks
f1e0088 Specify ilia stylesheet in all looks
```

### Changes in `regolith-i3xrocks-config`:
```
0b12b5a Use /proc/net/route for interface detection (#126)
```

### Changes in `regolith-ftue`:
```
e155b96 Remove unneeded sleep and launch shortcuts in all page mode
```

### Changes in `ilia`:
```
91e74dd Update docs
ab57873 Add page keybindings in all-page mode
050eac7 document page interface
e5e8a5a Add option -a for all page mode to address https://github.com/regolith-linux/ilia/issues/17.  Misc cleanup
cc80bc3 Fix icon in files page
a2f29f0 Change label to 'Commands' from 'Terminal'
3df2ca9 Add ability to launch host app to address https://github.com/regolith-linux/ilia/issues/18. Fix icon loading of notification view.
d451b54 Minor simplification to command launching
```

### Changes in `regolith-control-center`:
```
1be590f260444f80661db29d2566f072b5b69300 Remove unused power_saver_low_battery option in power panel
8310cf2bdd8463ce90b0e33f3f0aeffa8f76f3d4 Remove battery percentage toggle
0caac0d47 Update logo
```


## Changelog Delta from Regolith 1.6 to Regolith 2.0 Beta 1

### Changes in `regolith-lightdm-config`:
```
8a0d5c7 Switch font from jetbrainsmono to dejavusansmono
3eac0ff Update dependencies for Regolith 2
```

### Changes in `regolith-session`:
```
a939fad UI refresh logic moves to looks
cd72a50 gnome Xresource moves to looks
4f67559 Add Xsession file to add user desktop app dir to XDG_DATA_DIR to enable ilia to index desktop files there
e3c83d5 Remove unimplemented 'stage' command.  Misc cleanup
79cb6b1 Fix font spec from i3 to gnome
f776087 Only set gnome theme properties if defined by xres
870303f Set gnome interface font from look during refresh
9a9d98f Fix reference to term profile script and add soft dep to address problem in https://github.com/regolith-linux/regolith-look-extra/pull/1.
0ef3f55 Re-apply root Xresources after regolith Xres to ensure user values are set
a844012 Create new virtual package 'regolith-look-2' to specify a look. Version bump for cleaner upgrades from 1.x
d4bec14 Test that wallpaper file exists before setting
5c78271 Move Xresources.d from /etc to /usr/share/regolith due to apt policy issues with files in /etc.
2612147 Add breaking package from Regolith 1.6 due to file conflict
3a6b087 Fix sessoin init issue reported in https://github.com/regolith-linux/voulage/issues/2
c24ed0d Move to regolith2 config dir
db8953b Fix session init bug with invalid look path root
e89a150 Several bug fixes for Regolith 2 Xresource merging/loading
daf7368 Extract common session logic to lib file.  Fix regolith-look refresh and list
5f54d14 Add user defined post logoff script
614f6e3 Import regolith-look script from deprecated package regolith-styles.  Add build automation.
06f8c79 Fix path to renamed default look
eb85712 Depend on updated look package, add virtual look package.
```

### Changes in `regolith-i3-config`:
```
96c4762 Decompose navigation, workspace config, and base launchers from i3 config file.  Misc fixes
cf78c9d Add ability to specify custom stylesheet to ilia
8cad4c9 Remove unneeded dep from i3-gaps to i3-snapshot.  Cleanup
f8aec01 Remove session, resize, and user programs from hard dependencies in i3 config
0b7aa04 Extract i3 style config to optional partial
fa0e452 Allow resize mode keybinding to be overridden to address https://github.com/regolith-linux/regolith-desktop/issues/627
967fb85 remontoire and ilia keybindings conflict because they both specify help. Set 'conflicts' in debian accodringly.  Addresses https://github.com/regolith-linux/regolith-desktop/issues/628.
09eb0e2 Add label for gaps
72af37e Move i3 partials from /etc to /usr
e58e1d7 Move i3xrocks conf.d from /etc to /usr/share/regolith due to apt policy issues with files in /etc.
aa1ff5c Migrate from gnome-control-center to regolith-control-center
2b4e580 Fix ftue dependency
64fd601 Add breaking package from Regolith 1.6 due to file conflict
cdbb722 looks will host i3wm xres files, not here
fd0974e Move to regolith2 config dir
8e50f48 Add '--no-startup-id' when exec launchers
72a3636 Conditionally apply i3 xresources, allowing looks to override
b9072d5 Make rofi and ilia versions of rofication ui mutually exclusive
bf61321 Add seperate rofications for rofi and ilia front ends
1c23f65 Update alt variable to always occur on lhs due to i3 variable expansion constraints.  Related to https://github.com/regolith-linux/regolith-i3-gaps-config/pull/39
486151b Use ilia for look selector
118f557 use var for alt rather than hardcode. based from https://github.com/regolith-linux/regolith-i3-gaps-config/pull/39.
5880f98 Remove i3-next-workspace from common config
641f9d3 Switch to ilia for rofication UI
09765da Add all next-workspace bindings
6384067 add configurable workspace movement, from https://github.com/regolith-linux/regolith-i3-gaps-config/pull/37.
22815c6 Add keybinding viewer (ilia), replacing remontoire
73fda78 Add ilia i3 config
ac206fc Consolidate some partials to handle variable access limitations w/ i3 partials
d17ee44 build automation
2c5091f Update partials for i3 4.20
```

### Changes in `regolith-look-default`:
```
422a36b Specify custom stylesheet for ilia to fix Debian issue with incorrect font used in keybindings view
68041dc Wallpaper fixes to default loader.  Misc tweaks and cleanup
37ce868 Fix missing default_loader.sh in package
555f857 Fix bug when single background color is specified.
77d3d40 tweaks to font sizes.  Fix debian control
b2adcdb Refactor looks to be more consistent, own refresh logic, and better represent gtk theme elements
ed8cc44 Fix typeface for volume glyph
f6460f3 Create new virtual package 'regolith-look-2' to specify a look. Version bump for cleaner upgrades from 1.x
7381999 Alt char for temp
63f4fd1 Use alternate temp glphy to fix mystery ordering issue
0ccb167 Add missing icon theme dependency
5ef662a Move i3wm xres here. Fix non-active workspace color in bar
6d37a7f Fix glyphs in typeface not specified by font
1bf39d8 Update default look to use nerdfonts version of jetbrains mono.  update glyphs.
5b1973b Make looks non-exclusive
d778c94 Add rofication-ilia
0cc844b Switch from rofi to ilia
5cabe6e Fix gnome theme key
902a01b Add virtual package for a look
68d072b Change from minimal to default
```

### Changes in `regolith-look-extra`:
```
f1e0088 Specify ilia stylesheet in all looks
efe94ba Fatter borders in gruvbox for upcoming rounded corners
baaace4 Add hardy wallpaper to gruvbox. Publish gruvbox look
db7f92e Fix missed mapping from gtk monospace font to bar workspacce font. Reduce blackhole font size.  Assocated with https://github.com/regolith-linux/regolith-look-extra/issues/2
e572851 Add gruvbox look
7703e34 Tweaks polish and fixes to looks
066907f Sort glpyhs
67de712 Depend on default look for loader script
0ef9875 Merge branch 'main' of github.com:regolith-linux/regolith-look-extra into main
791321a Move blackhole to new Xres scheme
6944d7d Move nord to new Xres scheme
b8736b0 correct description and dependencies for blackhole
9da8e9c Package @ploum look blackhole.  Related https://github.com/regolith-linux/regolith-look-extra/pull/1
fdbbd77 Fix font on wifi and workspace glyph to address https://github.com/regolith-linux/regolith-desktop/issues/631
4420393 Fix font on volume indicator to address https://github.com/regolith-linux/regolith-desktop/issues/631
8c7656f Create new virtual package 'regolith-look-2' to specify a look. Version bump for cleaner upgrades from 1.x
3cddd15 More tweaks to term colors for nord
1f4bd7f Tweaks to term colors for nord
094fc90 visual refresh to nord look
9c956e1 Move nord to use nerdfonts version of sourcecodepro.  Misc fixes for Regolith 2 xres
965eac1 use original i3-wm Xres def for nord.  Disable gruvbox package.
4503a80 Fix def names for Regolith 2
0203187 Make looks non-exclusive
e2c9a3a Remove redundant look. Add nord adapted from 1.x
2e15aa8 Add docs, changelog
acb101f Initial commit. Ported from regolith-look-default
```

### Changes in `regolith-i3xrocks-config`:
```
7b2e609 Fix typo in media player script
ad20e67 i3xrocks Xresource moves into looks
9f73ef5 Specify font for both label and value after transition to nerdfonts
674e498 Fix celcius degree postfix on temp
87ab5d2 Migrate configuration from /etc to /usr/share due to problematic apt policy in /etc directories
885be5f Swap gnome-control-center for regolith-control-center
e8fd320 Add missing dep on rofication block
9170923 Fix ilia launch from notification
1b24e7b Deprecate python rofication status block for shell script. Launch ilia on click by default
```

### Changes in `regolith-desktop`:
```
2bf8c77 Push look dependencies down to regolith-look
f2f4435 Fix typo in packaging, use proper symver
e82a954 Revise package architecture.  Remove -small desktop package, go back to single base package.
12e39b6 Regolith 2 package model
adf1a53 Sync changelog to distro
a7c6444 Update versioning
```

### Changes in `regolith-ftue`:
```
6c5e76f Cleanup, better font default
a7bbfbb Add delay when loading ftue keybindings view in attempt to avoid race w/ i3 config loading
820a297 Use new conf directory for regolith2
1b542ed Add missing dependency
1135876 Switch to ilia over remontoire for keybinding viewer
0be6a48 Rename term profile script for consistency
eefcb55 Depend on uuid-runtime to fix debian ftue
```

### Changes in `regolith-compositor-xcompmgr`:
```
593c3e0 Move to regolith2 config dir
```

### Changes in `regolith-compositor-compton-glx`:
```
94cf43a Move to regolith2 config dir
6a234bf Add new name for zoom window, from https://github.com/regolith-linux/regolith-compositor-compton-glx/issues/3
454c5d6 Do not fade ilia
```

### Changes in `picom`:
```
e96a7b2 Package v9.0 release
```

### Changes in `regolith-compositor-picom-glx`:
```
49ddcf9 Move to regolith2 config dir
2ee5b27 Add new name for zoom window, from https://github.com/regolith-linux/regolith-compositor-compton-glx/issues/3
d27bbab Do not fade ilia
9fe5b87 Update to use @Filip-K picom config from https://github.com/regolith-linux/regolith-desktop/issues/578.
```

### Changes in `remontoire`:
```
f72a575 Cleanup source formatting
96090d5 Add scrollable parent to flowbox to allow vertical scrolling of content when exceeds bounds of flowbox
dc210de Don't draw scroll hints for horizontal scrolling
```

### Changes in `ilia`:
```
8443061 Tolerate empty css override. Add class name to keybindings view
b7ffa94 Numerous tweaks to help facility. Misc cleanup
6954aa1 Do not hardcode font in stylesheet
6f19ef1 Cleaner simpler approach to help content
761c5a8 Customize style via hardcoded CSS. Enable command line options for specifying css file or no styling.  Misc fixes.  Addresses https://github.com/regolith-linux/ilia/issues/12 and https://github.com/regolith-linux/ilia/issues/8
57b2037 Add ability to specify custom stylesheet
7ac9bda Add help feature, enabled by default
f64cc2e Better error message for invalid page
6967488 Fix array size for pages
66f6127 Update heuristic to grab resources at start time (experimental)
6130e48 Use XDG_DATA_DIRS to source desktop files to address https://github.com/regolith-linux/ilia/issues/11
83081f0 Load config partials for keybindings
e0a0c46 Update dev branches for automation
ec11206 Only grab keyboard and mouse
f5ac6f2 terminal page cleanup
0612e62 Add feature to resize dialog
499a084 Fix window launcher to close after action.  Misc cleanup
c4b0658 Add help flag to print all page types
9bc1b97 Add sleep after app launch to prevent intermittent failures
fa5163f Retry to grab inputs if fails
47dfe97 Use GTK POPUP window to avoid fullscreen window management of i3.  Relates to https://github.com/i3/i3/issues/4749
c1434d0 Fix textlist return selected value
cfa437d remove/disable debug output
a72ffc4 Add tracker dependency
8c040ca Tracker page functionally complete
5b8ff8f WIP: tracker partially working.
35637ec Update build files and add tracker dep
f67ede6 Window selection causes focus
00af209 Read windows from i3 and display
843bef4 Add partial parser for i3 tree struct
2adf24e Add working textlist page
97b1234 Fix message truncation bug
a979865 Do not handle any key events
4f3d557 Enable delete notification function
425c1b1 Allow pages to handle key events
e5af321 Icon and layout for notifications
f006839 parse i3 config file to extract default values for variables
ec55468 Keybindings viewable
961eda0 autovala updates
611fdb0 Add classes from remontoire to read i3 config
a76e39a WIP: icon/detail view.  waiting for model changes from rofication before proceeding
082f659 Can read/display notifications
1fc0c56 Add ability to specify primary page on command line
fec9a81 Use hashtable for app launch frequency lookups
3e93feb Prefer apps that start with query string.  Misc cleanup.
05c0989 Add UI polish
6450a6f WIP: terminal page
f0c2d4a Add spinner UI to clean up app launch lifecycle
5dc43f0 WIP, bug w/ intermittant launch failure
4e2db83 Extract desktop app loading into class.  Introduce general interface for pages.
b51bf40 Add app metadata files. Add build support for gsettings schema.
6517de8 Improve startup ux. Load apps from user dirs.  Close dialog before launching app.
35a1d08 Ignore whitespace in filter
111241d Support pg up/down. Configurable window and icon size.
8c78649 Add unused cmake files to appease autovala/meson
bbe2851 Selected apps rise to top of list
25a447d add gsettings schema
579c8d6 Can launch apps
b58252a Initial import from autovala creation
```

### Changes in `lago`:
```
960c8a0 Refactor rofication client to be sessionless.
6313183 More legible output from ls
31f495e Add empty files to appease autovala generated build files
b5bda93 Add usage to readme
1117188 Add readme
084b41d Add packaging
57939e8 clean up cli
43419b5 format source
03132a4 Add source format config
f85a450 fix src filename
bf452d8 Add rm and clear functions
fd4d765 CP: retrieve something from rofication dameon
```

### Changes in `regolith-rofication`:
```
a75b21a Fix ImportError exception for Python 3.10 or higher (#20)
7481617 Add command to delete list of notifications by id (#17)
6d4673b Merge pull request #16 from cfsmp3/master
cd1ecc8 Fix parsing hints
dab5dae Change path of queue storage file so it's not directly in the home dir
fb3d631 Merge pull request #15 from cfsmp3/master: Add timestamps to notifications
86285dc Undo unneeded changes in daemon
2be67af Made internal format a Python timestamp, convert to string in the UI
21b65d5 Fix conflicts
0ae8cca Read timestamp format from XResources
faeb17f Feature: store the notification state of 'icon' (#13)
3f14e05 Fix some typos in README (#14)
ffdc26e Add daemon
83f5b74 Make datetime format configurable via arg
1ed370a Add timestamp to notifications
e45205f Fix some typos in README (#12)
02329b5 Package last change
29ff56c keep ui open if there are still notifications left (#10)
```

### Changes in `regolith-control-center`:
```
6a39ca015 Add debug info for invalid panel id, related to https://github.com/regolith-linux/regolith-control-center/issues/2
cffbd8e1b Remove dim screen toggle
db3357d3b Remove background panel.  See https://github.com/regolith-linux/regolith-desktop/issues/630
87aab1ec1 For better i3 interop, all modal dialogs close if they lose focus, preventing poor UX from hidden window
a17b0366d Swap ubuntu logo for regolith
b60229237 Remove unneeded lock panel
a5e34f04c Remove unneeded notifications panel
2cfbbb664 Remove unused fractional scaling
558309824 Remove ubuntu and connectivity panels for debian support
001a1e711 Remove patch that fails to build on debian
0f95095a0 Remove gnome-language-selector support, ubuntu specific
158dcd919 Remove whoopsie support, ubuntu specific
bff27cc66 Update dbus service name to prevent collisiosn with gnome-control-center
2660efbf5 Rename desktop entry to avoid collision with gnome-control-center II
0c21d0d4b Rename desktop entry to avoid collision with gnome-control-center
be157d991 Remove upstream file
5d934dc83 Rename cc-remote-login-helper to prevent collision w/ gnome-control-center
17421fbbd Ignore files not installed in removed packages
44784f21d Remove unshipped files
a7a5a9fb7 Fix references to renamed schema
aa55ad004 Remove -data and -faces packages, use upstream
185af0ddd Rename gnome-power-manager icons to prevent collision w/ gnome-control-center
ffc1ed0f7 Update flatpack fields for regolith
a9859968f Change filenames to avoid collision w/ gnome-control-center
8bd341a1d Tweak desktop file to differentiate from gnome-control-center
695219bc7 Change filenames to avoid collision w/ gnome-control-center
024a179f3 Move gettext files due to collision with gnome-control-center
4a3db8b67 Add context to readme
8c1b2c9fd Remove failing dep
55c3a8d0f Remove failing or redundant patches
b160acbd3 Remove upstream readme
9b5bdb062 Regress polkit-gobject-1 dependency to match debian version in meson (experimental)
20bd1f69e Global rename from gnome-control-center to regolith-control-center
a86c8ef7a Fix of issue with adjusting output volume
```

### Changes in `regolith-default-settings`:
```
5aac9cb Remove texture from default wallpaper
```


### Changes in `nordic`:
```
3577075 Add built themes from release 2.1.0
dbd6e31 Remove source files
6f1aea6 Add debian packaging
```

### Changes in `fonts-source-code-pro-ttf`:
```
2fd5dfd Initial import from speedricer, debian cleanup
```

### Changes in `fonts-nerd-fonts`:
```
96d4985 Remove noto as assets are too large for github policy
a511489 add source file
a63a1f9 remove script from debian branch
7842791 Commit from automation
```
