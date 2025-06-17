---
title: EpsteinLib
layout: hextra-home
toc: false
date: 2025-02-10T16:27:07+02:00
math: true
---




<div class="hx-container">


{{< hextra/hero-container
  image="/images/zeta.svg"
  imageTitle="EpsteinLib lattice sum example"
  imageWidth="400"
>}}


{{< hextra/hero-badge
link="https://github.com/epsteinlib/epsteinlib/blob/main/CHANGELOG.md"
_comment="Can we Link html page here ?"
>}}
  <div class="hx-w-2 hx-h-2 hx-rounded-full hx-bg-primary-400"></div>
  <span>Latest version: v1.0</span>
  {{< icon name="arrow-circle-right" attributes="height=14" >}}
{{< /hextra/hero-badge >}}

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-headline >}}
  EpsteinLib
{{< /hextra/hero-headline >}}
</div>

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-subtitle >}}
High-performance C Library for the computation of the Epstein zeta function
{{< /hextra/hero-subtitle >}}
</div>

<div class="hx-mt-6 hx-mb-6">
{{< hextra/hero-button text="Get Started" link="docs/" >}}
</div>
{{< /hextra/hero-container >}}


<div class="hx-mt-6"></div>
<div class="hx-mt-6"></div>
<div class="hx-mt-6"></div>



{{< hextra/feature-grid cols="2" >}}

{{< hextra/feature-card
  title="Accurate"
  subtitle="Full Accuracy over the whole parameter range. <br> <img src=\"images/error.svg\" style=\"max-height:500px; width:auto; display:block; margin: 20px auto 0;\" />"
  style="background: radial-gradient(ellipse at 50% 80%,rgba(194,97,254,0.15),hsla(0,0%,100%,0));"
>}}

{{< hextra/feature-card
  title="Fast"
  subtitle="3D Lattice sums in less than $0.01$ milliseconds. <br> <img src=\"images/timing.svg\" style=\"max-height:500px; width:auto; display:block; margin: 20px auto 0;\" />"
  style="background: radial-gradient(ellipse at 50% 80%,rgba(142,53,74,0.15),hsla(0,0%,100%,0));"
>}}

{{< /hextra/feature-card >}}

<div class="hx-mt-6"></div>

{{< hextra/feature-grid cols="1" >}}
  {{< hextra/feature-card
    title="Accessible"
    link="https://pypi.org/project/epsteinlib/"
    link="/docs"
    subtitle="Get the  Python wrapper with `pip install epsteinlib` or follow the <span class=\"regolith-links\">installation</span> instructions."
style="background: radial-gradient(ellipse at 50% 80%,rgba(221,210,59,0.15),hsla(0,0%,100%,0));"
  >}}
{{< /hextra/feature-grid >}}

<div class="hx-mt-6 hx-mb-6"></div>
<div class="hx-mt-6 hx-mb-6"></div>
{{< hextra/hero-section >}}
  Examples
{{< /hextra/hero-section >}}


{{< columns cols="2" >}}

  {{< column
      title="Trivial zeros of the Epstein zeta function"
      border="true"
      image="/images/zeros.svg"
  >}}
Works like this.
{{< /column >}}

  {{< column
      title="Anomalous Quantum Dispersion Relation"
      border="true"
      image="/images/dispersion.svg"
  >}}
Works like this.
{{< /column >}}

{{< column
    title="Singular Euler--Maclaurin expansion"
    border="true"
    image="/images/sem.png"
>}}
Works like this.
{{< /column >}}

{{< column
    title="Casimir Effect"
    border="true"
    image="/images/casimir.png"
>}}
Works like this.
{{< /column >}}


</div>


{{< /columns >}}
