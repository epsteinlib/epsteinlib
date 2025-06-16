---
authors:
  - name: Khosrow Moossavi
    link: https://github.com/khos2ow
    image: https://avatars.githubusercontent.com/u/372575?v=4
date: 2025-02-01
linktitle: Archive URL change for upcoming v3.3 release
title: Archive URL change for upcoming v3.3 release
---

For the Regolith v3.3 release, the host name of the URL from which `apt` pulls packages as well as the structure of package archive repositories have changed.  The new hostname is `archive.regolith-desktop.com`.  In order to update to Regolith v3.3 and get fixes and upgrades going forward kindly update your apt configuration.

The old `apt` URL:

```text {filename="/etc/apt/sources.list.d/regolith.list"}
deb [arch=amd64] https://regolith-desktop.org/release-current-ubuntu-noble-amd64 noble main
```

The new `apt` URL:

```text {filename="/etc/apt/sources.list.d/regolith.list"}
deb [arch=amd64] https://archive.regolith-desktop.com/ubuntu/stable noble main
          |                                           |      |      |     |
          |                                           |      |      |     * Repo Component
          |                                           |      |      * Distro Version/Codename
          |                                           |      * Repo Suite
          |                                           * Distro Name
          * Architecture
```

{{< callout type="info" >}}
All other aspects of the package configuration remain the same including the key used to verify packages.
{{< /callout >}}

### Architecture

| Name | Description |
|------|:------------|
| `amd64` | The 64-bit x86 architecture |
| `arm64` | The 64-bit Arm architecture |

### Distro Names

| Name | Description |
|------|:------------|
| `ubuntu` | The Ubuntu Linux Distribution |
| `debian` | The Debian Linux Distribution |

### Repo Suites and Components

| Regolith Stage | Repo Suite | Repo Component | Description |
|----------------|:-----------|:--------------:|:------------|
| Experimental | `experimental` | `main` | Stage for testing experiments |
| Unstable | `unstable` | `main` | Stage for early development testing |
| Testing | `testing` | `main` | Stage for late development testing |
| Release X.Y | `stable` | `vX.Y` | Stage for specific official release |
| Current Release | `stable` | `main` | Stage for latest release (floating from release to release) |

### Distro Version/Codename

These labels are determined by their respective upstream communities.  Examples are `noble`, `bookworm`, `focal`.

### Examples

| Description | APT Line |
|:------------|:---------|
| The 3.2 release of Regolith on Ubuntu Jammy for `amd64` | `deb [arch=amd64] https://archive.regolith-desktop.com/ubuntu/stable jammy v3.2` |
| The latest release of Regolith on Debian Bookworm for `amd64` | `deb [arch=amd64] https://archive.regolith-desktop.com/debian/stable bookworm main` |
| Testing repo for Regolith on Debian Testing for `arm64` | `deb [arch=arm64] https://archive.regolith-desktop.com/debian/testing testing main` |
