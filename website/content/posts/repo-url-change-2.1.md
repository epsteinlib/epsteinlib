---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2022-08-04
linktitle: Package URL change for 2.1
title: Package URL change for 2.1
---

## Summary

For the Regolith 2.1 release, the host name of the URL from which `apt` pulls packages has changed.  The new hostname is `regolith-desktop.org`.  In order to update to Regolith 2.1 and get fixes and upgrades going forward kindly update your apt configuration.

The old `apt` URL:

```text {filename="/etc/apt/sources.list.d/regolith.list"}
deb [arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] https://regolith-desktop.io/release-ubuntu-jammy-amd64 testing main
                                                                                                     ^
                                                                                                     |
                                                                                                     (this changes to '.org')
```

The new `apt` URL:

```text {filename="/etc/apt/sources.list.d/regolith.list"}
deb [arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] https://regolith-desktop.org/release-ubuntu-jammy-amd64 testing main
```

All other aspects of the package configuration remain the same including the key used to verify packages.

## Additional Info

The Regolith 2.1 release reintroduced the ISO form of distribution.  In order to generate the ISO with the production `apt` URL the packages needed to be published first.  This essentially required a way of publishing the release before we wish users to get access to the packages.  In order to achieve this the idea of using a URL and switching the host seemed like a good idea.  Problems occurred when the realities of TLS cert management became apparent.  When the domain was remapped to the new host, the TLS certificates were not.  Later, when the problem became known, the hosting provider gave advice on how to resolve the problem.  However, the original source domain name (`regolith-desktop.io`) somehow is not able to take a new cert registration from the new host.  Some remedies have been suggested by users, but due to the nature of the production host we opted to leave the system as-is and migrate to the new domain to prevent further churn for users.

When the 2.1 release was published, the URL was updated to `regolith-desktop.org`.  Users still using `regolith-desktop.io` will most likely be stuck on the 2.0 release unless the TLS certificate issue is ignored or in some way resolved.
