---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2023-11-26
linktitle: Status Updates Nov 2023
title: Status Updates Nov 2023
---

## Regolith 3.1 & Ubuntu Mantic Support

Due to an [issue](https://github.com/regolith-linux/regolith-desktop/issues/936) associated with the Desktop Portal system and the GNOME implementation, that results in a severe usability regression, Ubuntu Mantic support has been delayed until it can be resolved.  This change came into Regolith with the inclusion of GNOME 45 into Ubuntu Mantic.  It looks like we're getting some help from upstream, but it may be some time until we're able to support GNOME 45-based distros if the solution ends up being that we need to implement our own portal backend to coexist with GNOME's backend.

Additionally, Ubuntu 24.04 will be an LTS release and the people behind Ubuntu are putting extra work into it to ensure it's a great one.  We may end up forgoing Mantic support altogether if a short term solution does not present itself.

## Ubuntu Summit 2023

Soumya and I had the opportunity to travel to Riga, Latvia earlier this month to attend the [Ubuntu Summit](https://ubuntu.com/blog/ubuntu-summit-2023), graciously hosted by Canonical.  There were a lot of great sessions and it was fantastic to reconnect with some familiar faces and collaborators from the previous year. So many brilliant and friendly people in one place!

![](/images/posts/P_20231103_185135~2.jpg)

In particular [PapyElGringo](https://github.com/PapyElGringo) announced a new project [Veshell](https://github.com/free-explorers/veshell) as a successor to the [Material Shell](https://github.com/material-shell/material-shell) project.  I enjoyed learning about the spacial approach to app and workspace management, and look forward to seeing how the new project progresses.

 Soumya and I attended the Ubuntu Flavor sync session and were given several helpful contacts to work with as the Regolith project gets a closer association to upstream projects Debian and Ubuntu.  Stay tuned for future news on this front.

## Chat System Migration: Slack -> Matrix

Long-time Regolith contributor Leo had setup a Slack instance in the early days of the Regolith project.  It was a great help and it was an effective tool for collaboration.  Over time though, I and others in the project had a desire to move to an open platform.  In part this is simply due to the risky nature of relying on "free" versions of commercial, proprietary platforms.  Partly due to the Ubuntu's adoption of the [Matrix protocol](https://matrix.org/) for some community collaboration, I found it to be functional and easy to use.  A [Matrix "Space" for the Regolith project](https://matrix.to/#/#regolith-desktop:matrix.org) has been setup, any everyone is welcome to join.  At some point in the next month or so we'll retire the Slack server entirely.
