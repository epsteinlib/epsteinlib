---
title: "Install"
linkTitle: "Install"
weight: 1
description: >
  Install Regolith
prev: /docs/using-regolith
---

## Upgrading to 3.x

If you're upgrading from an earlier version of Regolith, you may be interested in reading [the release notes](/docs/reference/Releases/regolith-3.2-release-notes) and if coming from 2.x, the [migration guide](/docs/reference/Releases/regolith-3.0-release-notes/#migration-guide).

## Supported Operating Systems

### Ubuntu

{{< tabs items="Ubuntu 24.04,Ubuntu 22.04" >}}
{{< tab >}}
Regolith can be installed as system packages.  This makes updating and removing easier and more consistent.  The following steps describe how
to configure your system to read packages from the Regolith package repository and install the desktop package.

1. Register the Regolith public key to your local `apt`:

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. Add the repository URL to your local `apt`:

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/ubuntu/stable noble v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. Update `apt` and install Regolith

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. System Restart

The login manager will need to be restarted for the new desktop session to be recognized. The easiest way of restarting it is to reboot your system.

{{< callout type="info" >}}
Replace `amd64` with `arm64` in the two places in the above line to install on ARM-based systems.
{{< /callout >}}
{{< /tab >}}

{{< tab >}}
Regolith can be installed as system packages.  This makes updating and removing easier and more consistent.  The following steps (also available as a [downloadable script](/install-release-ubuntu-22.04-amd64.txt)) describe how
to configure your system to read packages from the Regolith package repository and install the desktop package.

1. Register the Regolith public key to your local `apt`:

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. Add the repository URL to your local `apt`:

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/ubuntu/stable jammy v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. Update `apt` and install Regolith

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. System Restart

The login manager will need to be restarted for the new desktop session to be recognized. The easiest way of restarting it is to reboot your system.

{{< callout type="info" >}}
Replace `amd64` with `arm64` in the two places in the above line to install on ARM-based systems.
{{< /callout >}}
{{< /tab >}}
{{< /tabs >}}

### Debian

{{< tabs items="Debian Testing,Debian Bookworm" >}}
{{< tab >}}
For Debian users that use the "Testing" release, Regolith can be installed, however there is no "release" version due to the nature of Debian Testing.  In order to install Regolith into a Debian Testing instance, the `apt` line contains `testing` instead of a `release` string.  For example:

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/debian/testing testing main" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

Otherwise the installation steps are the same as an official Debian release.  Note that broken dependencies or components will occur from time to time.  If stability is important, consider installing into an officially released version of Debian or Ubuntu.
{{< /tab >}}

{{< tab >}}
Regolith can be installed as system packages.  This makes updating and removing easier and more consistent.  The following steps (also available as a [downloadable script](/install-release-debian-12-amd64.txt)) describe how
to configure your system to read packages from the Regolith package repository and install the desktop package.

1. Register the Regolith public key to your local `apt`:

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. Add the repository URL to your local `apt`:

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/debian/stable bookworm v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. Update `apt` and install Regolith

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. System Restart

The login manager will need to be restarted for the new desktop session to be recognized. The easiest way of restarting it is to reboot your system.

{{< callout type="info" >}}
Replace `amd64` with `arm64` in the two places in the above line to install on ARM-based systems.
{{< /callout >}}

{{< callout type="info" >}}
The `regolith-desktop` package installs a reasonable minimal configuration that should work on most hardware. To add additional package for a richer and more functional environment read [Recommended Packages for New Users](/docs/using-regolith/configuration/#recommended-packages-for-new-users).
{{< /callout >}}

{{< callout type="info" >}}
The `regolith-compositor-picom-glx` compositor should work on most computers. If you experience driver or visual issues, [try another compositor](docs/howtos/customize-compositor).
{{< /callout >}}
{{< /tab >}}
{{< /tabs >}}

## Packages Specified at Install Time

Starting with Regolith 3.0 it is now possible to select from multiple sessions (X11: `regolith-session-flashback`, Wayland: `regolith-session-sway`). Due to this, when installing Regolith 3.2 also specify one or both sessions that you wish to use. Additionally, it's now possible to specify the Regolith Look at install time, which saves a setup of having to configure it later.  Here is the recommended base install for X11 on Debian-based systems:

```bash
sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
#                ^-- base package ^-- session                ^ -- look 
```

* For the session, you may chose one or both of: `regolith-session-flashback` (X11), `regolith-session-sway` (Wayland)
* For the look, you may chose *one* from any available look package:
```text
regolith-look-ayu-dark
regolith-look-ayu-mirage
regolith-look-ayu
regolith-look-blackhole
regolith-look-default-loader
regolith-look-default
regolith-look-dracula
regolith-look-gruvbox
regolith-look-i3-default
regolith-look-lascaille
regolith-look-nevil
regolith-look-nord
regolith-look-solarized-dark
```

An an example and alternative to the recommendation above, here is the `apt` line to install the Sway session with the Nord Look:

```bash
sudo apt install regolith-desktop regolith-session-sway regolith-look-nord
```

## Anatomy of the Regolith `apt` Config Line

```text
deb [arch=amd64] https://archive.regolith-desktop.com/ubuntu/unstable noble main
          |                                           |      |        |     |
          |                                           |      |        |     * Repo Component
          |                                           |      |        * Distro Version/Codename
          |                                           |      * Repo Suite
          |                                           * Distro Name
          * Architecture
```

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
|:---------------|------------|:--------------:|:------------|
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

## Release Change Policy

From Regolith 3.0 to 3.2 (inclusive), all releases will use unique names in the package repo URL.  This means that users will be in full control of when they wish to upgrade to a new release.  Users wishing to have their package manager always install the latest version, a special stage called `release-current` is provided.

From Regolith 3.3 onward, all releases will be published in the common archive repository separated by a high level distro (e.g. `debian`, `ubuntu`) folder and corresponding codename subfolders. The structure is split furthermore into different suites and components. As such, users wishing to have their package manager always install the latest version, they should use `main` component of `stable` suite. Otherwise a fixed version component (for example `v3.3`, `v3.2`, etc) can be used alongside `stable` component.

Note that `experimental`, `unstable`, and `testing` suites only have `main` component.
