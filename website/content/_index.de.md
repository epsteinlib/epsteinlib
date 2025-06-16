---
title: "Die Regolith Desktop Umgebung"
layout: "hextra-home"
---

{{< hextra/hero-container
  image="/images/v-tour/regolith-empty.png"
  imageTitle="Die Regolith Desktop Umgebung"
  imageWidth="853"
>}}
{{< hextra/hero-badge link="/docs/reference/releases/regolith-3.2-release-notes" >}}
  <div class="hx-w-2 hx-h-2 hx-rounded-full hx-bg-primary-400"></div>
  <span>Latest version: v3.2</span>
  {{< icon name="arrow-circle-right" attributes="height=14" >}}
{{< /hextra/hero-badge >}}

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-headline >}}
  Die Regolith Desktop Umgebung
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
    title="Produktivität steht an erster Stelle"
    subtitle="Regolith verwendet i3: ein beliebter, schneller und konfigurierbarer Kachelfenster-Manager, der sich hervorragend für schnelle tastaturgesteuerte Arbeitsabläufe eignet. Regolith integriert i3 mit anderen Desktop-Komponenten wie z.B. `i3bar`, `rofication`, `gnome-flashback`, und `ilia`, um eine vollständige Desktop-Oberfläche zu bieten."
  >}}

  {{< hextra/feature-card
    title="GNOME System Management"
    subtitle="Durch die Bereitstellung einer Gnome-Sitzung mit `gnome-flashback`, umgeht Regolith einen Großteil der Komplexität einer typischen `gnome-shell`-basierten Umgebung und bietet gleichzeitig einfache und konsistente Systemverwaltungsmöglichkeiten."
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
    Bei der ersten Anmeldung wird ein Overlay (jederzeit aktivierbar mit
    {{< keys "super,?" >}}) die wichtigsten Tastenkombinationen, die in i3-wm
    verwendet werden.
  {{< /column >}}

  {{< column
      title="Terminal"
      border="true"
      image="/images/v-tour/regolith-floating-terminal.png"
  >}}
    Für diejenigen, die im Terminal arbeiten, genügt die Tastenkombination
    {{< keys "super,enter" >}} um loslegen zu können.
  {{< /column >}}

  {{< column
      title="App Launcher"
      border="true"
      image="/images/v-tour/regolith-ilia-apps.png"
  >}}
    Der globaler App-Launcher ist jederzeit verfügbar, um Ihre Programme mit
    {{< keys "super,space" >}} zu starten.
  {{< /column >}}

  {{< column
      title="Manage Layout"
      border="true"
      image="/images/v-tour/regolith-desktop-terminals.png"
  >}}
    Sie brauchen mehr Terminals? Erstellen Sie Layouts im Handumdrehen, indem
    Sie zwischen horizontalem und vertikalem Modus umschalten mit
    {{< keys "super,backspace" >}} gefolgt von {{< keys "super,enter" >}}.
    Navigieren Sie zu Fenstern mit {{< keys "super,h,j,k,l" >}}.
  {{< /column >}}

  {{< column
    title="Floating Window"
    border="true"
    image="/images/v-tour/regolith-floating-windows.png"
  >}}
    Wechseln Sie den Modus für schwebende Fenster mit {{< keys "super,f" >}}.
    Ändern der Fenstergröße geht mit {{< keys "super,r" >}} und bewegen Sie sie
    mit der Maus, indem Sie {{< keys "super" >}} gedrückt halten.
  {{< /column >}}

  {{< column
      title="Gnome System Managament"
      border="true"
      image="/images/v-tour/regolith-gnome-flashback.png"
  >}}
    Gnome Flashback bietet eine konsistente und einfache Systemverwaltung. Passen
    Sie Ihre Benutzeroberfläche an, binden Sie Ihre USB-Laufwerke automatisch ein,
    verbinden Sie sich mit drahtlosen Netzwerken. Starten Sie die Systemsteuerung
    mit {{< keys "super,c" >}}.
  {{< /column >}}

  {{< column
      title="Regolith Look"
      border="true"
      image="/images/v-tour/regolith-screenshot-look-selector.png"
  >}}
    Wechseln Sie einfach zu <a href="https://ethanschoonover.com/solarized" class="regolith-links">Solarized</a>
    oder anderen looks mit dem <code>regolith-look</code> Befehl. Da
    <b>looks</b> den Paketmanager verwenden, speichern Sie nur die Themenressourcen,
    die Sie gerade verwenden.
  {{< /column >}}

  {{< column
      title="Find and Navigate"
      border="true"
      image="/images/v-tour/regolith-ilia-windows.png"
  >}}
    Haben Sie viel um die Ohren? Finden Sie schnell das gesuchte Fenster über
    {{< keys "super,ctrl,space" >}} oder navigieren Sie über Arbeitsbereiche mit
    {{< keys "super,[0 - 19]" >}}.
  {{< /column >}}

  {{< column
      title="Manage Gaps"
      border="true"
      image="/images/v-tour/regolith-many-windows.png"
  >}}
    Verschwenden Sie keinen Platz für eine überflüssige Benutzeroberfläche und
    nutzen Sie jeden Pixel, ohne zum Mikro-Manager Ihrer Fensterlayouts zu werden.
  {{< /column >}}

  {{< column
      title="Notifications"
      border="true"
      image="/images/v-tour/regolith-ilia-notifications.png"
  >}}
    Desktop-Benachrichtigungen kämpfen nicht um Ihre Aufmerksamkeit, sondern
    können über einen Dialog auf dem Bildschirm verwaltet werden, mit den
    Tasten {{< keys "super,n" >}}.
  {{< /column >}}
{{< /columns >}}
