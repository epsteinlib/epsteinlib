---
title: "Regolith Desktop Environment"
layout: "hextra-home"
---

{{< hextra/hero-container
  image="/images/v-tour/regolith-empty.png"
  imageTitle="Regolith Desktop Environment"
  imageWidth="853"
>}}
{{< hextra/hero-badge link="/docs/reference/releases/regolith-3.2-release-notes" >}}
  <div class="hx-w-2 hx-h-2 hx-rounded-full hx-bg-primary-400"></div>
  <span>Latest version: v3.2</span>
  {{< icon name="arrow-circle-right" attributes="height=14" >}}
{{< /hextra/hero-badge >}}

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-headline >}}
  Regolith Desktop Environment
{{< /hextra/hero-headline >}}
</div>

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-subtitle >}}
  Productivity focused Ubuntu derivative and desktop environment
{{< /hextra/hero-subtitle >}}
</div>

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-button text="Get Started" link="docs/" >}}
{{< hextra/hero-badge style="padding: 13px 12px !important; font-size: .875rem !important;" link="https://opencollective.com/regolith/donate" >}}
  <span>Support Regolith <img class="not-prose" style="display: inline; height: 22px;" src='https://badgen.net/opencollective/backers/regolith'/></span>
  {{< icon name="arrow-circle-right" attributes="height=14" >}}
{{< /hextra/hero-badge >}}
</div>
{{< /hextra/hero-container >}}

<div class="hx-mt-6"></div>
<div class="hx-mt-6"></div>
<div class="hx-mt-6"></div>

{{< hextra/feature-grid cols="2" >}}
  {{< hextra/feature-card
    title="Productivity First"
    subtitle="Regolith runs i3 and Sway: popular, fast, and configurable tiling window managers that are great for productive, keyboard-driven workflows. Regolith integrates i3 and Sway with other desktop components such as `i3bar`, `rofication`, `i3status-rs`, and `ilia` to provide a complete desktop experience."
  >}}

  {{< hextra/feature-card
    title="GNOME System Management"
    subtitle="By integrating light-weight window managers with `gnome-session`, Regolith sidesteps much of the complexity of a typical `gnome-shell` based environment while providing simple and consistent system management capabilities."
  >}}
{{< /hextra/feature-grid >}}

<div class="hx-mt-6"></div>

{{< hextra/feature-grid cols="1" >}}
  {{< hextra/feature-card
    title="Install Regolith"
    link="/docs/using-regolith/install/"
    subtitle="Follow Regolith Desktop Installation instructions <span class=\"regolith-links\">here</span>."
  >}}
{{< /hextra/feature-grid >}}

<div class="hx-mt-6 hx-mb-6"></div>
<div class="hx-mt-6 hx-mb-6"></div>
{{< hextra/hero-section >}}
  Virtual Tour
{{< /hextra/hero-section >}}

{{< columns cols="2" >}}
  {{< column
      title="First Login"
      border="true"
      image="/images/v-tour/regolith-ilia-keybindings.png"
  >}}
    Upon first login, an overlay (toggled anytime via {{< keys "super,?" >}})
    presents the most important keybindings used with i3-wm.
  {{< /column >}}

  {{< column
      title="Terminal"
      border="true"
      image="/images/v-tour/regolith-floating-terminal.png"
  >}}
    For those that work in the terminal, pressing {{< keys "super,enter" >}}
    is all it takes to get to business.
  {{< /column >}}

  {{< column
      title="App Launcher"
      border="true"
      image="/images/v-tour/regolith-ilia-apps.png"
  >}}
    A single global app launcher is instantly available from anywhere to run
    your programs via {{< keys "super,space" >}}.
  {{< /column >}}

  {{< column
      title="Manage Layout"
      border="true"
      image="/images/v-tour/regolith-desktop-terminals.png"
  >}}
    Need more terminals? Create layouts on the fly by toggling between horizontal
    and vertical modes with {{< keys "super,backspace" >}}, followed by
    {{< keys "super,enter" >}}. Navigate to windows positionally with
    {{< keys "super,h,j,k,l" >}}.
  {{< /column >}}

  {{< column
    title="Floating Window"
    border="true"
    image="/images/v-tour/regolith-floating-windows.png"
  >}}
    Toggle floating window mode with {{< keys "super,f" >}}. Resize windows with
    {{< keys "super,r" >}} and move them around with the mouse by pressing
    {{< keys "super" >}}.
  {{< /column >}}

  {{< column
      title="Gnome System Managament"
      border="true"
      image="/images/v-tour/regolith-gnome-flashback.png"
  >}}
    Gnome Flashback provides consistent and simple system management. Tweak your
    UI, auto mount your USB drives, connect to wireless networks. Launch the
    control panel with {{< keys "super,c" >}}.
  {{< /column >}}

  {{< column
      title="Regolith Look"
      border="true"
      image="/images/v-tour/regolith-screenshot-look-selector.png"
  >}}
    Easily switch to <a href="https://ethanschoonover.com/solarized" class="regolith-links">Solarized</a>
    or other looks with the <code>regolith-look</code> command. Because **looks**
    utilize the package manager, you only store theme resources that you're using.
  {{< /column >}}

  {{< column
      title="Find and Navigate"
      border="true"
      image="/images/v-tour/regolith-ilia-windows.png"
  >}}
    Got a lot going on? Quickly find the window you're looking for via
    {{< keys "super,ctrl,space" >}} or navigate over workspaces with
    {{< keys "super,[0 - 19]" >}}.
  {{< /column >}}

  {{< column
      title="Manage Gaps"
      border="true"
      image="/images/v-tour/regolith-many-windows.png"
  >}}
    Waste no space on frivolous UI and take advantage of every pixel without
    micro-managing your window layouts.
  {{< /column >}}

  {{< column
      title="Notifications"
      border="true"
      image="/images/v-tour/regolith-ilia-notifications.png"
  >}}
    Desktop notifications do not compete for your attention, but rather can be
    managed via an on-screen dialog by pressing {{< keys "super,n" >}}.
  {{< /column >}}
{{< /columns >}}
