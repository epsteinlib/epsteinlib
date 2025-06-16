---
title: "Build System"
linkTitle: "Build System"
weight: 8
description: >
  General information regarding how Regolith packages and ISOs are built
next: /docs/reference/releases
---

## Overview 

This page describes how packages are generate from source code, and how those packages combine to actualize the Regolith desktop environment.  This document covers Regolith only and does not intend to describe upstream dependencies, general package manager concepts, or other aspects of build systems.

## Terms

* `Distro` - Short for distribution, in the context of this document, may be considered a variation of a Linux-based OS operating system.  A product, project, and brand.
* `Codename` - A well defined and announced snapshot of a set of packages that are consumed as a whole.  May be considered a software release with a specific name such as `36` or `Jammy`.
* `Architecture` - In the context of this document, denotes a CPU and base system hardware.  Regolith currently targets the `amd64` and `arm64` architectures.
* `Stage` - Represents a snapshot of a set of packages at a modeled level of maturity.  In Regolith the following stages from least-stable to most-stable are modeled:
  1. `experimental` - For packages and package changes that may be reverted, removed, or change significantly over time.  Serious bugs are common and expected here.
  2. `unstable` - For packages changes that are planned to reach a release but have yet to undergo much testing.  Serious bugs may occur from time to time.  Regolith core developers may run their testing and development from this stage but users are discouraged from consuming packages from `unstable`.
  3. `testing` - For packages undergoing testing by Regolith developers and early adopters.  Bugs occur but serious bugs should not often pop up in this stage.  Regolith users that wish to get access to fixes and new features earlier may consume from this stage with the understanding that bugs are to be expected (, and reporting them is appreciated!)
  4. `release` - For packages that have been determined to be suitable for general use.  While minor bugs may be present in packages at this stage, the bugs present should not interfere with general use of the environments for users on supported distros and releases.
* `Suite`  - Represents the published packages repository. It has a one-to-one corresponding relationship to `Stage`. In other words `Stage` has the point of view of software development and `Suite` has the point of view of published archive repositories.
  1. `experimental` suite - corresponds to `experimental` stage
  2. `unstable` suite - corresponds to `unstable` stage
  3. `testing` suite - corresponds to `testing` stage
  4. `stable` suite - corresponds to `release-X_Y` and `release-current` stages
* `Component`  - Each published repository is divided into several components. They correspond to Regolith released versions. In Regolith these are `main` and correposnding versions (for example `v3.2`, `v3.1`, etc)
* `Target` ~ Represents a specific sequence of `distro`, `codename`, `stage`, and `architecture` for which a set of packages exist.
* `Manifest` ~ A file that specifies a set of the properties name, origin, and snapshot for a list of packages.  A manifest is a file that represents a snapshot of a `target` in time.
* `Package Manager` ~ A set of formats and programs that allow for packages to be installed in a distro.  The package manager may very by distro.
* `Package Repository` ~ A well defined source of binary packages that may be installed into an OS that supports a given package manager.


## Source Repositories

All Regolith packages are generated from git-versioned source repositories.  These repositories house the source code for a given component.  Mostly, each component is general for any given `distro` and `codename`.  Source repositories are modeled in the Regolith build system via `manifest` files.  

## Manifests

The Regolith build system reads and writes text files known as manifests.  The files are initially read as the source of truth for what a given `target` package repository should contain.  If the build system determines if there are updates to the manifest present in the manifest's source repositories, then any updates are built into package updates and published to the target repo along with the updated manifest.  In the case that no updates are found, the build simply completes without any changes.  As such, a given target's manifest can be used to determine what packages, versions, and source branches/tags are present.

### Manifest Schema

Manifests are simply text files with one package per line, and each element in on a line separated by a space (ASCII 32).

```text
<PACKAGE NAME> <PACKAGE REPO URL> <PACKAGE REPO BRANCH/TAG> <PACKAGE REPO COMMIT HASH>
```

## Package Model

The Regolith build system utilizes sets of JSON files to define what packages are specified for a given `target`.  These JSON files are set in a specific structure to facilitate overriding higher-scoped contexts for lower-scoped contexts.  For example, at the top level is the `release` directory.  It specifies all packages for that release.  Within a `release`, directories for specific `distro`s exist.  Package models in these directories may override or nullify packages specified from it's parent.  In this way, a top-level general package model can be tweaked for a given `stage`/`distro`/`codename`/`architecture`.  To further illustrate this point here is the structure of the `stage` directory in the [`voulage` repo](https://github.com/regolith-linux/voulage) (which is the name of the Regolith 2.x package builder) (NOTE: some files were removed for brevity):

```text
stage/
├── package-model.json                <--- Root package model
├── release
│   ├── debian
│   │   └── bullseye
│   │       └── package-model.json
│   ├── package-model.json            <--- Package model for the release stage
│   └── ubuntu
│       └── jammy
│           └── package-model.json
├── testing
│   ├── debian
│   │       └── package-model.json
│   ├── package-model.json
│   └── ubuntu
│       └── jammy
│           └── package-model.json
└── unstable
    ├── debian
    │   ├── bullseye
    │   │   └── package-model.json
    │   └── testing
    │       └── package-model.json    <--- Customizations for Debian Testing in unstable stage
    ├── package-model.json
    └── ubuntu
        └── jammy
            └── package-model.json

```

## Packages

As mentioned above, the build system generates packages for a given target from it's manifest.  The build system is designed to be agnostic to any particular package format.  As such, package versions are not directly modeled in manifests.  The mechanism that determines if a source change results in a version bump is delegated to package-manager specific extensions.  After a package has been identified to require building, the build system delegates to the package-specific build system to generate the target packages, such as debian source and binary packages.  If the package update builds successfully, another package-manager extension is called to inject the newly built package into the distro-dependant package repository.

## Package Build System

The Regolith build system consists mainly of bash [shell scripts](https://github.com/regolith-linux/voulage/.github/scripts).  These scripts are run from [GitHub workflows](https://github.com/regolith-linux/voulage/.github/workflows) but are designed to be runnable directly from the appropriate local environment to facilitate easy testing and troubleshooting.  The following is a high level summary of what the package builder does once invoked: 

1. Set up source state such as the target string, which is generated from `<stage>-<distro>-<codename>-<architecture>`
2. Get the current manifest from the Voulage repository (inside `manifests/` folder)
3. Generate the merged package model starting from the root and working into each segment in the `target` string
4. With the merged model, iterate over each package, generating a manifest entry by resolving the commit hash from the specified ref
5. After all manifest entries are generated, compare the generated manifest with the earlier downloaded copy
6. If no changes found, exit because no packages have changes on the specified refs
7. Otherwise, package changes have been detected.  Download the entire target repository
8. For each difference found, delegate to distro-specific build package build and repo publish commands
9. Publish the new packages to the target repository
9. Commit and push the updated manifest back to the Voulage repository

## ISO Build System

The ISO builder is based on the open source project [live-custom-ubuntu-from-scratch](https://github.com/mvallim/live-custom-ubuntu-from-scratch/).  The ISO builder uses a set of shell scripts to generate an Ubuntu ISO image for both a live environment and an installer.  Regolith customizes this project in [this repo](https://github.com/regolith-linux/regolith-ubuntu-iso-builder) with some specific packages to [install and remove](https://github.com/regolith-linux/regolith-ubuntu-iso-builder/blob/main/scripts/regolith-2_1-jammy-config.sh).  A github workflow is used to generate the ISOs and upload them to GitHub for hosting.

## Cookbook

### How to add or update a package to the build system

1. Determine which `stage`s, `distro`s, `codename`s, and `architecture`s the change applies to (ex: `unstable-ubuntu-jammy-amd64`)
2. Edit the package model in the `stage` directory in the [Voulage](https://github.com/regolith-linux/voulage)
3. Add or edit the package stanza JSON for that package ([example](https://github.com/regolith-linux/voulage/commit/cd4a3ebed9c501610810e7d4c79c4811c992c426))
4. Commit the change
5. Execute the [build workflow](https://github.com/regolith-linux/voulage/actions)
6. Once the build completes, check the `releases` in Voulage to verify that the generated packages were produced and published as expected

### How to determine what commit a package is built with in a given target package repo

1. Get the manifest contained in the Voulage repository (located at `manifests/<DISTRO>/<CODENAME>/<SUITE>-<COMPONENT>/<ARCH>/manifest.txt`)
2. Find the line containing the package name
3. Note the branch and ref in the manifest, this is the repo commit from which the package was built

### How to extend the build system for a new package format

1. For an example read the [Debian v4 build workflow](https://github.com/regolith-linux/voulage/blob/main/.github/workflows/build-deb-v4.yml)
2. Refer to the [Debian extension script](https://github.com/regolith-linux/voulage/blob/main/.github/scripts/ext-debian.sh) that the workflow uses
3. Create an extension script that maps from the git source of each package to the target distro's packaging format
