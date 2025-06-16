---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2022-02-05
linktitle: Announcing Regolith 2
title: Announcing Regolith 2
---

## Introduction

Regolith 2 is the next major release of the Regolith Desktop environment. Here is a summary of the big changes in the release:

* Use `i3` or `i3-gaps` interchangeably
* Simplify the loading of Xresources, discovery, and management of looks
* Support for i3 config partials, allowing for safer and richer extension capabilities
* Provide a package topography necessary for future changes including alternative sessions and backends
* Debian Bullseye support, along with Ubuntu Focal and Jammy
* Move to new custom desktop executor `ilia` from `rofi` for app and command launching, notifications, and more

Regolith 1.6 is the end of the line for the existing Regolith desktop. The new major version 2 provides for a simplification of Regolith as a whole, with some new features as well.  This upgrade will not be backwards compatible with 1.x series (but as always, users may choose to specify their own i3 config). However Regolith 2 is not radically different than the previous version. This major version bump is to allow for users to upgrade on their own schedules.  As such, the 1.x package repository (PPA) will not be used for Regolith 2 packages.

## Release Status

We're currently in pre-alpha phase of the development process. The environment is generally stable for day to day work, however because we're reserving the right to make big changes still, using it in a production context is not recommended at this time.

## Shipped Changes

* Convert i3 config to [partials](https://github.com/i3/i3/issues/4192)
* `conf.d` style Xresource discovery
* Migration from `rofi` to [`ilia`](https://github.com/regolith-linux/ilia) as the desktop executor
* Migration to [`nerdfonts`](https://www.nerdfonts.com/) for looks to simplify and enhance font icons
* A [new build system](https://github.com/regolith-linux/voulage) decoupled from `launchpad.net` facilitating easier and more consistent distro development
* Port of looks `lascaille` and `nord` from the 1.x looks

## Changes in Development

* New website (this one!)
* [Polish](https://github.com/regolith-linux/ilia/issues) for `ilia`
* Bugs in wallpaper handling (gnome/dotfiles)
* Better HiDPI support out of the box
* Better first time user experience
* Porting most existing looks to new Xresource schema
* Rebranding from "Regolith Linux" to "Regolith Desktop Environment"

## Expected Release Date

Due to the volunteer nature of the Regolith project, it's simply not possible to a promise ship date.  However, we would like to have something usable by the time of the next LTS version of Ubuntu is released, in April 2022.

## Expected Upgrade Experience

The Regolith 2 desktop environment will require the removal of the existing Regolith 1.x PPA and the addition of a new Debian repository.  `regolith-desktop` will continue to be the top-level package that provides the complete desktop environment.  As before, i3 config files can be entirely controlled by users by staging a copy in their XDG user configuration directory, however with `i3` `4.20+`, a new option of simply providing any additional configuration via config partials is also supported.

## Learn More

[See here for details](http://regolith-desktop.com/#get-regolith-20) on how to install the current pre-release version of Regolith 2, and [here for information on how you can provide feedback or contribute](http://regolith-desktop.com/#get-regolith-20).
