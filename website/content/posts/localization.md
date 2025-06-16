---
authors:
  - name: Ken Gilmer
    link: https://github.com/kgilmer
    image: https://avatars.githubusercontent.com/u/49683
date: 2022-06-22
linktitle: Looking for Translators
title: Looking for Translators
---

With the 2.0 release out the door, we're looking for help to translate the index page and the "using regolith" section into languages other than English. In the past we've used a translation system, but it turned out to be too complex for many people to understand. For the [Regolith 2.0 site](https://github.com/regolith-linux/regolith-desktop-website) we'll rely on [Hugo's multilingual support](https://gohugo.io/content-management/multilingual/).

For a given language translation the following files must be translated:

```text
content/_index.md
content/docs/using-regolith/basics.md
content/docs/using-regolith/configuration.md
content/docs/using-regolith/first-launch.md
```

To provide a translation, simply copy the English file and add the [two-letter language code](https://lingohub.com/academy/best-practices/iso-639-1-list) to the end of the filename.  Once the translation is complete, create a pull request and we'll merge your translation.  This allows other people to read about Regolith with your language translation.

As an example, if a Japanese translation was provided, the following new files would be added to the PR:

```text
content/_index.ja.md
content/docs/using-regolith/basics.ja.md
content/docs/using-regolith/configuration.ja.md
content/docs/using-regolith/first-launch.ja.md
```

There is some additional Hugo configuration before the translation will be available, but we can take care of that.  Feel free to open an issue or reach out on [GitHub Discussions](https://github.com/orgs/regolith-linux/discussions) with any questions, and thanks in advance!

## Translation Workflow

Here are some step-by-step instructions for creating a translation:

1. Create a fork of the [GitHub repo that contains the Regolith Desktop website](https://github.com/regolith-linux/regolith-desktop-website).  ([learn more](https://docs.github.com/en/get-started/quickstart/fork-a-repo))
2. Clone the fork to your local computer.
3. Determine the [two-letter designation](https://lingohub.com/academy/best-practices/iso-639-1-list) for your language.
4. Copy the following files to new files with the two-letter language designation (see example above).
5. Translate the content to your language.
6. Commit and push the changes to your fork of the website repository.
7. Create and submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).

By creating the pull request, Regolith project members will be notified to review and push your changes.
