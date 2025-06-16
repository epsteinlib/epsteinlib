# This is the https://regolith-desktop.com website

This site uses the [Hugo static site generator](https://gohugo.io/) and the [Hextra theme](https://themes.gohugo.io/themes/hextra/).

## Local Development

Please ensure you have the following installed:

- Hugo **extended** **>= 0.140.0** (otherwise PostCSS won't work)
- npm
- git

### Checkout and View Local

Clone this repo, and the Book theme as a submodule:

```console
$ git clone git@github.com:regolith-linux/regolith-desktop-website.git regolith-desktop-website
$ cd regolith-desktop-website
```

Install dependencies and launch Hugo server:

```console
$ hugo
$ hugo server
```

...or, with Docker:

```console
docker compose up
```

### Make changes

Edit content as needed. If `hugo server` is running then updates are automatically reloaded after committing changes to disk.

### Update website on GitHub Pages

GitHub Actions on the repository is configured to deploy any and all changes to the `gh-pages` branch automatically, no further action needed. The `CNAME` ([https://regolith-desktop.com](https://regolith-desktop.com)) is configured automatically, no need to add it.

If you believe that your changes ought to be reviewed by a second party, just open a pull request and have somebody else jump in and cross-check your edits.

### Shortcodes

Additional shortcodes are:

- `keys`: To insert a recognizable sequence of keypresses, e.g. `{{< keys "super,shift,q" >}}`. The block also accepts an optional delimiter, should you need `,` as a keypress, e.g. `{{< keys "super|," "|" >}}`
- `img`: To insert an image, including automated, cross-device sizing, a lightbox and an alt text, e.g. `{{< img "images/my-image.jpg" "Some alt text" >}}`
- `supported_version`: To insert the Ubuntu versions, delimited by commas, which are currently supported by Regolith Linux, e.g. `{{< supported_versions >}}` (they are taken from the `config.toml` `params` section). You can either choose to display the version _numbers_ or the version _names_, with _numbers_ being the default. To display the names use `{{< supported_version type="names" >}}`

### Conventions

#### Use Markdown wherever possible

Although Hugo accepts common HTML elements pretty much anywhere the expressed goal is to rely on Markdown and its syntax for pretty much anything _within the confines of the `content` directory itself_. If you need to style something differently, try using [a shortcode](https://gohugo.io/templates/shortcode-templates/) first, before actually styling the page in "pure" HTML. Some pages are violating this concern, mostly because the underlying theme does as well, but for the sake of making it easier on folks to edit and especially translate the content as well as keeping a consistent style across languages they should be free of HTML.

#### Use `prettier` for formatting files (especially Markdown)

[`prettier`](https://www.npmjs.com/package/prettier) is a versatile formatter for JavaScript/TypeScript and everything related to HTML as well as JSON/YAML and TOML. Please use it to format all the content in the repository that `prettier` is capable of processing. Manual formatting will show up as an error in the future and PRs changing formatting that's different from `prettier` output will not be accepted.

#### Use `console` for highlighting shell session code

If you want to show a command that users can run, you are encouraged to use `bash` for code highlighting. But if the command you are showing has an output you need to use `console`.

This block uses `bash` for highlighting:

```bash
sudo apt update
sudo apt install regolith-desktop
```

But this uses `console` (note that this command is also showing and output):

```console
$ xrandr --listmonitors
Monitors: 1
 0: +*eDP-1 1920/309x1080/173+0+0  eDP-1
```

Using `bash` obscures the flow of information and things you type into a session isn't really shell code all the time but rather a random set of commands. There's [an open feature request](https://github.com/alecthomas/chroma/issues/137) with the highlighting engine ("Chroma") Hugo is using to accept `console` as a lexer. Once that's accomplished all the relevant highlighting parts will look different. For now, they'll just receive the default styling.


#### Call it `Xresources` (with backticks)

No `Xresource` or Xresource. Use the plural (because that's what we read and what the definition of the mechanism is) and surround the word itself with backticks to make it stand out a little.

#### Put images into `static/images` and use the `img` shortcode

We don't want duplicate content anywhere, especially not images that are the same for every single translation. Therefor we put images into the `static/images` folder and use the `img` shortcode to refer to them. You'll get a lightbox and automated, screen-aware (i.e. "mobile friendly") sizing on top. The `img` shortcode has the following syntax:

```go-html-template
{{< img "<path-to-image-in-asset-folder>" "alt-text-for-image" >}}
```

Example:

```go-html-template
{{< img "images/my-image.jpg" "Some image" >}}
```

#### Use `ref` or `relref` for linking between content

Content changes, page names change. Using a hardcoded path in Markdown links will break sooner or later. Hugo [has two built-in functions](https://gohugo.io/content-management/cross-references/), `ref` (preferred) or `relref`, to make sure links are viable and linking to the right content. Use them whenever you're linking to "internal" content.

**Note**: Hugo will fail to render the page if a link is "dead" and you will have to fix it before you can continue. We consider this a feature!

#### Use the Markdown footnote syntax to render footnotes

Hugo allows for specifying footnotes in [the regular Markdown syntax](https://michelf.ca/projects/php-markdown/extra/#footnotes), e.g. `[^1]` for the note in the text and `[^1]: Some text` (mind the colon) wherever you want to place the footnote text. It has the added advantage that the footnote also links back to the text.
