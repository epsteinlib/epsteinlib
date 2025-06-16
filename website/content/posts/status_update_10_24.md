---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2024-10-29
linktitle: Status Updates October 2024
title: Status Updates October 2024
---

## Ubuntu 24.10 (Oracular Oriole) Support

Canonical completed their 20th Ubuntu release with Oracular this month.  A beta of Regolith 3.2 is available for testing, see [this discussion thread](https://github.com/orgs/regolith-linux/discussions/1073) for details.  As this is not a new Regolith release, but rather the addition of a new Ubuntu version on Regolith 3.2, the [3.2 release notes](/docs/reference/Releases/regolith-3.2-release-notes/) continue to be the latest.

## Ubuntu Summit 2024

<img src="/images/posts/P_20241027_083737.jpg" alt="Near Huis ten Bois at Ubuntu Summit 24" title="Near Huis ten Bois at Ubuntu Summit 24" width="320" style="float: right; margin: 35px;"/>

I submitted a [lightning talk](https://events.canonical.com/event/51/contributions/554/) that was accepted at [Ubuntu Summit 2024](https://events.canonical.com/event/51/) in [The Hague](https://en.wikipedia.org/wiki/The_Hague), Netherlands.  I was able to attend in person and while present at the conference, had a wonderful time touring around a few of the many many parks by bike. The talk itself was intended to be a simple project status update. However during the conference a number of [inspiring](https://events.canonical.com/event/51/contributions/509/) [talks](https://events.canonical.com/event/51/contributions/522/) motivated me to try for more.  In particular, it struck me during a demo of the new Cosmic Desktop, that rarely if ever did I hear anyone discussing *why* i3-like windowing environments are attractive to some users. So, I attempted to briefly explain my take during my brief talk:

> i3-like environments (as typically configured, of course you can do whatever you like) bias towards application usage over discovery.  Concretely, i3 expects a bit more of the user, to know what they want to do.  In return, the system throws less distraction and minutia at the user.  This, at lease in my experience, leads to longer periods of deep concentration on the problem I'm working on. In summary, the user is better able to put their attention where it matters: the application.

Due to the short period of the talk, I didn't have much time to get my point across, but hopefully the information was interesting to some attendees. I'm grateful for the event organizer's continued support of the Regolith Desktop project.

## Regolith 3.3

Work continues on the next Regolith release.  Highlights include:

* Simplified Deb repository layout
* `ilia` module refactor
* `armesto` integration
