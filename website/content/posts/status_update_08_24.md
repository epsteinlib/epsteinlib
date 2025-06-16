---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2024-08-02
linktitle: Status Updates August 2024
title: Status Updates August 2024
---

## Regolith 3.2 

We're currently on beta 4 of Regolith 3.2.  It's been a longer release cycle than hoped, as many users are wanting to run Regolith on the latest Ubuntu LTS: Noble.  In this release, the items that come to mind that required more time than expected were the systemd init refactoring and integration of GNOME 46.  Soumya initiated the effort to clean up and make our init code generally consistent with systemd, cleaning up some rather old ad-hoc shell script init code.  This refactoring covers both Sway and i3 window managers.  I expect this effort to provide a more robust and consistent start up process that in turn will allow users to more easily configure and reason about the execution of their desktop setup.

We will require at least one more beta to clean up some issues encountered between X11 and Wayland sessions, as well as take in some great fixes provided by community members.   But, it feels like we are almost there.

## DebConf 24 

<img src="/images/posts/PXL_20240729_102150414.jpg" alt="Soumya and Ken at DebConf 24" title="Soumya and Ken at DebConf 24" width="320" style="float: right; margin: 35px;"/>

Due to the generous support of our financial contributors, Soumya and I were able to meet up at [DebConf 24](https://debconf24.debconf.org/) this year in Busan, Korea.  This is our second in-person conference meet up; last year Canonical graciously sponsored both of us to attend the Ubuntu Summit in Riga, Latvia.  That event was incredibly productive and so we were both keen on further opportunities to collaborate in person and with the larger Debian world.

The travel costs associated with our attendance at the conference was paid from our [Open Collective account](https://opencollective.com/regolith/donate).  While we do not spend a lot of effort to advertise or solicit financial contributions, happily there are enough people donating to the project to allow us to travel from time to time when advantageous to the project and it fits into our personal schedules. There are also a few financial contributors over the years that have contributed rather generous sums to the project.  While I don't feel it appropriate to call them out by name here, we are very happy to have received these funds to support the project.  

We aim to invest all of the financial contributions back into the project such that users will receive value back in terms of better quality software, broader support, and healthy working relationships with related projects.

## What's Coming Next

### Debian Project Integration

Soumya has [begun the effort](https://packages.debian.org/sid/trawld) to integrate the Wayland session into the Debian package repository.  Our goal is to integrate the entire wayland session into Debian.  Once there, we plan to then provide a "spin" or "flavor" to downstream Ubuntu.  I hope that this will eventually result in an officially supported flavor that will enable users to install Regolith from an ISO as has been done in past LTS releases.  However due to time constraints, I will not take on providing an ISO for the 24.04 release.

### Simplified Debian Repositories

Khosrow has been working for a number of months on an [updated Debian apt repo setup](https://github.com/regolith-linux/voulage/pull/38).  Once this is in place, the Regolith deb repos will be far more consistent with other projects and will be easier for users to configure and reason about.

### Armesto

Work on a rust-based replacement to the [rofication daemon](https://github.com/regolith-linux/armesto) is expected to land in the next release, in addition to the Wayland session getting a more consistent waybar configuration.  This will also set the stage to enable some user-requested capabilities of notifications, such as a stronger visual indication of new notifications.

### Ilia

More work is planned to go into [`ilia`](https://github.com/regolith-linux/ilia) to modularize pages such that it is usable by any general desktop environment.  The Regolith-specific logic will move to plugins provided by distinct binary packages.  This will allow easier integration into Debian and provide another launcher option to the free desktop community.

Additionally I will continue to work on load performance improvements for users that experience lag in the dialog loading of applications.

## Deprecations

As mentioned in the 3.2 release notes, we will no longer produce packages for Debian Bullseye or any Ubuntu releases other than Jammy and Noble.  

We will not provide an ISO of the 24.04 release.

## Feedback

[Feedback](https://github.com/orgs/regolith-linux/discussions) about this status update or the Regolith project in general is always welcome.
