---
title: 'Change background'
linkTitle: 'Change background'
description: >
    Change the background
---

## Wallpaper

Most Regolith Looks provide a default desktop color or image. The wallpaper can be set by specifying the path to the wallpaper image via the Xresources key `regolith.wallpaper.file` or in the Regolith Control Center which can be opened by running `regolith-control-center` in the terminal.

Example, assuming the file `/usr/share/backgrounds/hardy_wallpaper_uhd.png` is present:

```bash
echo "regolith.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

In the above command the property mentioned is set to the `/usr/share/backgrounds/hardy_wallpaper_uhd.png` change this path to any other that needs to be used.

To change the image treatment (zoom, scale, etc.) of the wallpaper image, set `regolith.wallpaper.options`. Options are: `none`, `wallpaper`, `centered`, `scaled`, `stretched`, `zoom`, `spanned`

Example:

```bash
echo "regolith.wallpaper.options: zoom" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

To specify a color rather than an image use the Xresource key `regolith.wallpaper.color.primary`:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

To specify a secondary color and gradient:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.secondary: green" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.shading.type: vertical" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Lockscreen Wallpaper

The lockscreen wallpaper can be managed in the same way as the desktop wallpaper. The keys are the same as desktop wallpaper, but with `lockscreen` prefixed to `wallpaper`, like so: `regolith.lockscreen.wallpaper...`. For example, to specify a lockscreen image:

```bash
echo "regolith.lockscreen.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### Disable Wallpaper handling

If you wish to manage wallpaper externally to Regolith, simply specify empty values for wallpaper image and color:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: " >> ~/.config/regolith3/Xresources
```
