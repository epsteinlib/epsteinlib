---
title: "L'environnement Regolith Desktop"
layout: "hextra-home"
---

{{< hextra/hero-container
  image="/images/v-tour/regolith-empty.png"
  imageTitle="L'environnement Regolith Desktop"
  imageWidth="853"
>}}
{{< hextra/hero-badge link="/docs/reference/releases/regolith-3.2-release-notes" >}}
  <div class="hx-w-2 hx-h-2 hx-rounded-full hx-bg-primary-400"></div>
  <span>Latest version: v3.2</span>
  {{< icon name="arrow-circle-right" attributes="height=14" >}}
{{< /hextra/hero-badge >}}

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-headline >}}
  L'environnement Regolith Desktop
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
    title="La productivité d'abord"
    subtitle="Regolith utilise i3: un gestionnaire de fenêtre en tuiles qui est     populaire, rapide, configurable et très pratique pour tout piloter avec son clavier. Regolith intègre i3 avec d'autres composants comme `i3bar`, `rofication`, `gnome-flashback`, et `ilia` afin de fournir une interface de bureau complète."
  >}}

  {{< hextra/feature-card
    title="Système de gestion GNOME"
    subtitle="En fournissant une session gnome avec `gnome-flashback`, Regolith s'affranchi d'une bonne partie de la complexité des environnements typiques basés sur `gnome-shell` en fournissant un système de gestion simple et consistant."
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
    Lors de la première connexion, un dialogue (activé à n'importe quand avec
    {{< keys "super,?" >}}) présente les raccourcis clavier les plus important
    utilisés avec i3-wm.
  {{< /column >}}

  {{< column
      title="Terminal"
      border="true"
      image="/images/v-tour/regolith-floating-terminal.png"
  >}}
    Pour celles et ceux qui travaillent dans le terminal, pressez
    {{< keys "super,enter" >}} suffit à vous satisfaire.
  {{< /column >}}

  {{< column
      title="App Launcher"
      border="true"
      image="/images/v-tour/regolith-ilia-apps.png"
  >}}
    Un lanceur d'application global est immédiatement accessible depuis partout
    pour lancer vos programmes avec {{< keys "super,space" >}}.
  {{< /column >}}

  {{< column
      title="Manage Layout"
      border="true"
      image="/images/v-tour/regolith-desktop-terminals.png"
  >}}
    Besoin de plus de terminaux ? Créez un agencement à la volée en basculant
    entre mode horizontal et verical avec {{< keys "super,backspace" >}} suivi
    de {{< keys "super,enter" >}}. Naviguez vers les fenêtres selon leur position
    avec {{< keys "super,h,j,k,l" >}}.
  {{< /column >}}

  {{< column
    title="Floating Window"
    border="true"
    image="/images/v-tour/regolith-floating-windows.png"
  >}}
    Activez le mode fenêtre flottante avec {{< keys "super,f" >}}. Redimmensionnez
    les fenêtres avec {{< keys "super,r" >}} et déplacez-les avec la souris en
    pressant {{< keys "super" >}}.
  {{< /column >}}

  {{< column
      title="Gnome System Managament"
      border="true"
      image="/images/v-tour/regolith-gnome-flashback.png"
  >}}
    Gnome Flashback fourni un système de gestion simple et consistant. Personnalisez
    votre UI, montez automatiquement vos disques USB, connectez vous à des réseaux
    sans fils. Vous pouvez lancer le panneau de contrôle avec {{< keys "super,c" >}}.
  {{< /column >}}

  {{< column
      title="Regolith Look"
      border="true"
      image="/images/v-tour/regolith-screenshot-look-selector.png"
  >}}
    Changer facilement de thème pour <a href="https://ethanschoonover.com/solarized" class="regolith-links">Solarized</a>
    ou d'autres thèmes avec la commande <code>regolith-look</code>. Comme les
    **thèmes** utilisent le gestionnaire de paquet, vous stockez seulement les
    ressources des thèmes que vous utilisez.
  {{< /column >}}

  {{< column
      title="Find and Navigate"
      border="true"
      image="/images/v-tour/regolith-ilia-windows.png"
  >}}
    Trop de choses sur le feu? Trouvez rapidement la fenêtre que vous cherchez
    avec {{< keys "super,ctrl,space" >}} ou naviguer entre les espaces de travail
    avec {{< keys "super,[0 - 19]" >}}.
  {{< /column >}}

  {{< column
      title="Manage Gaps"
      border="true"
      image="/images/v-tour/regolith-many-windows.png"
  >}}
    Ne perdez pas de place avec des interfaces frivoles et profitez de chaque
    pixel sans avoir à faire du micro-management de votre agencement de fenêtres.
  {{< /column >}}

  {{< column
      title="Notifications"
      border="true"
      image="/images/v-tour/regolith-ilia-notifications.png"
  >}}
    Les notifications de bureau ne se battent pas pour votre attention mais
    peuvent être gérées depuis un dialogue en pressant {{< keys "super,n" >}}.
  {{< /column >}}
{{< /columns >}}
