---
title: "Konfiguration"
linkTitle: "Konfiguration"
weight: 4
description: >
 Nehmen Sie Änderungen am Aussehen und Verhalten von Regolith vor.
next: /docs/howtos
---

## Empfohlene Pakete für neue Benutzer

### Looks

Installieren Sie alle offiziellen Regolith-Looks:

```bash
sudo apt install regolith-look-*
```

### Status-Indikatoren

Das folgende Snippet installiert Statusindikatoren-Pakete, die für neue Benutzer empfohlen werden:

```bash
sudo apt install i3xrocks-focused-window-name i3xrocks-rofication i3xrocks-info i3xrocks-app-launcher i3xrocks-memory
```

Um einen dieser Indikatoren später zu entfernen, deinstallieren Sie ihn einfach von Ihrem System. Beispiel:

```bash
sudo apt remove i3xrocks-info
```

#### Laptop-Akkustatus

Zeigt den aktuellen Batteriestatus eines Laptops an:

```bash
sudo apt install i3xrocks-battery
```

## Hintergrundbild

Die meisten Regolith Looks bieten eine Standard-Desktop-Farbe oder ein Bild. Das Hintergrundbild kann durch Angabe des
Pfades zum Hintergrundbild über den Xresources-Parameter `regolith.wallpaper.file` festgelegt werden.

Beispiel: Angenommen, die Datei `/usr/share/backgrounds/hardy_wallpaper_uhd.png` ist vorhanden:

```bash
echo "regolith.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh 
```

Um die Darstellung (Zoom, Skalierung, etc.) des Hintergrundbildes zu ändern, setzen Sie `regolith.wallpaper.options`.
Optionen sind: `none`, `wallpaper`, `centered`, `scaled`, `stretched`, `zoom`, `spanned`

Beispiel:

```bash
echo "regolith.wallpaper.options: zoom" >> ~/.config/regolith3/Xresources
regolith-look refresh 
```

Um eine Farbe anstelle eines Bildes anzugeben, verwenden Sie den Xresource-Parameter `regolith.wallpaper.color.primary`:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
regolith-look refresh 
```

Um eine sekundäre Farbe und einen Farbverlauf zu definieren:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.secondary: green" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.shading.type: vertical" >> ~/.config/regolith3/Xresources
regolith-look refresh 
```

### Sperrbildschirm-Hintergrundbild

Das Hintergrundbild des Sperrbildschirms kann auf die gleiche Weise wie das Desktop-Hintergrundbild verwaltet werden.
Die Parameter sind die gleichen wie beim Desktop-Hintergrundbild, aber mit dem Präfix `lockscreen` vor `wallpaper`, etwa
so: `regolith.lockscreen.wallpaper...`. Um zum Beispiel ein Bild für den Sperrbildschirm anzugeben:

```bash
echo "regolith.lockscreen.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh 
```

### Deaktivieren der Wallpaper-Verwaltung

Wenn Sie das Hintergrundbild außerhalb von Regolith verwalten möchten, geben Sie einfach leere Werte für das
Hintergrundbild und die Farbe an:

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: " >> ~/.config/regolith3/Xresources
```

## Statusleisten-Anzeigen

Statusanzeigen wie CPU-Auslastung, Datum und Uhrzeit, Benachrichtigungen, Wetter und andere Systeminformationen können
durch die Installation von Paketen hinzugefügt oder entfernt werden. Um zum Beispiel eine Statusanzeige für die Batterie
Ihres Laptops anzuzeigen, führen Sie einfach `sudo apt install i3xrocks-battery` aus und aktualisieren Sie dann die
Sitzung mit `regolith-look refresh`. Um herauszufinden, welche Indikatoren verfügbar sind, führen
Sie `apt search ^i3xrocks-` aus oder suchen Sie nach `i3xrocks-` in Ihrer bevorzugten Paketmanager-GUI, wie zum
Beispiel [Synaptic](https://help.ubuntu.com/community/SynapticHowto). Es gibt [weitere Dokumentation]({{< ref "/docs/howtos/add-remove-blocklets.md" >}}) für die Konfiguration von Statusanzeigen, auch "blocklets" genannt.

## Looks

Farben, Hintergrundbilder, Fenster- und Leistenlayouts und andere visuelle Faktoren werden in Regolith gebündelt und
als "Looks" bezeichnet. Looks bieten eine einfache Möglichkeit, das gesamte Aussehen des Desktops zu ändern. Wie die
Statusleistenindikatoren werden auch die Looks in Paketen zusammengefasst und wie jedes andere Softwarepaket installiert
und entfernt. Look-Pakete haben das folgende Namensformat: `regolith-look-<name>`. Um nach verfügbaren Looks zu suchen,
kann `apt` oder ein GUI-Paketmanager verwendet werden. Zur Auswahl eines installierten Looks wird mit der
Tastenkombination {{< keys "super,alt,l" >}} ein Dialog angezeigt:

![](/images/v-tour/regolith-ilia-look-selector.png)

Alternativ kann der Terminalbefehl `regolith-look` verwendet werden, um das Aussehen zu ändern und die aktive Sitzung
mit dem ausgewählten Aussehen zu aktualisieren. Hier ist ein Beispiel, das zum `gruvbox`-Look wechselt:

```console
$ apt search ^regolith-look-
[...]
regolith-look-gruvbox/unknown,now 0.4.6-1regolith amd64
[...]
$ sudo apt install regolith-look-gruvbox
regolith-look set gruvbox
regolith-look refresh
```

Um einfach alle verfügbaren Looks zu installieren:

```bash
sudo apt install regolith-look-*
```

## i3 Features

Seit Regolith 2.0 werden viele Aspekte der i3-Konfiguration über den Paketmanager verwaltet. Mit i3 `4.20` wurde die
Möglichkeit eingeführt, Include-Dateien als Teil der Konfiguration anzugeben. Regolith stellt seine gesamte
i3-Konfiguration über diese Dateien bereit, die auch "config partials" genannt werden. Diese Pakete installieren ihre
i3-Konfigurationspartials in `/usr/share/regolith/i3/config.d`. Durch das Installieren und Entfernen von Paketen kann
die i3-Konfiguration an die eigenen Vorlieben angepasst werden, während gleichzeitig die Upstream-Änderungen für Aspekte
der Konfiguration, die sich nicht ändern müssen, unverändert bleiben können. Wenn das Paket `regolith-desktop`
installiert wird, werden diese Konfigurationselemente standardmäßig auch als Soft-Dependencies installiert:

| Paket                        | Funktion          |
|------------------------------|-------------------|
| regolith-i3-base-launchers   | Terminal und Browser starten | 
| regolith-i3-default-style    | Fensterverhalten |
| regolith-i3-navigation       | Navigationstastenbelegung |
| regolith-i3-networkmanager   | Netzwerk- und WiFi-Funktionen |
| regolith-i3-next-workspace   | Arbeitsbereich-Dienstprogramm |
| regolith-i3-resize           | Tastaturbelegung Größenänderung |
| regolith-i3-session          | Tastaturbelegung Session |
| regolith-i3-unclutter        | Mauszeiger ausblenden, wenn er nicht benutzt wird |
| regolith-i3-user-programs    | Optionales Starten von Benutzerprogrammen, die in Xresources angegeben sind |
| regolith-i3-workspace-config | Tastaturbelegung des Arbeitsbereichs |

Soft-Abhängigkeiten können entfernt werden, ohne dass Pakete entfernt werden, die davon abhängen. Dies bedeutet, dass
jedes der aufgelisteten Pakete entfernt werden kann. Benutzer können Änderungen an der Konfiguration auf solidere Weise
vornehmen, indem sie einen Teil der Standardkonfiguration durch ihre eigene Version ersetzen. Dies kann erreicht werden,
indem man den Teil, der angepasst werden soll, nach `~/.config/regolith3/i3/config.d/` kopiert und die Originaldatei
in `/usr/share/regolith/i3/config.d` über apt entfernt. Zum Beispiel, um die Tastaturbelegung des Arbeitsbereichs
anzupassen:

```bash
mkdir -p ~/.config/regolith3/i3/config.d
cp /usr/share/regolith/i3/config.d/40_workspace-config ~/.config/regolith3/i3/config.d/
vim ~/.config/regolith3/i3/config.d/40_workspace-config # gewünschte Änderungen vornehmen
sudo apt remove regolith-i3-workspace-config # Entfernt die Standardversion
# Starte i3 neu oder logge dich wieder ein
```

Auf diese Weise können die unveränderten Teilkonfigurationen weiterhin Updates und Fehlerbehebungen erhalten, ohne die
benutzerspezifische Konfiguration zu beeinflussen.

### Alle i3-Konfigurationspakete

Im Folgenden finden Sie eine Liste aller i3-Konfigurationspakete, die in Regolith 2.0 verfügbar sind:

| Paket                        | Funktion                                                                  |
|------------------------------|---------------------------------------------------------------------------|
| regolith-i3-base-launchers   | Terminal und Browser starten                                              | 
| regolith-i3-compositor       | Compositor-Integration (delegiert an `regolith-compositor-<variant>`)     |
| regolith-i3-default-style    | Fensterverhalten                                                          |
| regolith-i3-ftue             | Hilfe für neue User                                                       |
| regolith-i3-gaps-partial     | `i3-gaps` Konfiguration (wenn installiert, ersetzt i3-gaps i3)            |
| regolith-i3-gnome            | Regolith-Einstellungen                                                    |
| regolith-i3-i3xrocks         | Konfiguration der i3-Balkenstatusanzeigen                                 |
| regolith-i3-ilia             | `ilia` Einbindung                                                         |
| regolith-i3-navigation       | Tastaturbelegung für die Navigation                                       |
| regolith-i3-networkmanager   | "netzwerk-manager"-Einbindung                                             |
| regolith-i3-next-workspace   | Funktion zum automatischen Erstellen des nächsten Arbeitsbereichs         |
| regolith-i3-remontoire       | Integration der Ansicht "Legacy keybindings"                              |
| regolith-i3-resize           | Tastaturbelegung für Größenänderung                                       |
| regolith-i3-rofi             | Integration des alten Desktop-Startprogramms                              |
| regolith-i3-rofication       | Integration der alten Benachrichtigungsanzeige                            |
| regolith-i3-rofication-ilia  | Integration der Benachrichtigungsanzeige                                  |
| regolith-i3-session          | Tastaturbelegung der Session                                              |
| regolith-i3-snapshot         | i3-Arbeitsbereichspersistenz                                              |
| regolith-i3-unclutter        | Ausblenden der Maus bei Nichtgebrauch                                     |
| regolith-i3-user-programs    | Optionaler Start von Benutzerprogrammen, die in Xresources angegeben sind |
| regolith-i3-workspace-config | Arbeitsbereich-Tastaturbelegung                                           |

## Tastaturbelegungen

{{< callout type="warning" >}}
Regolith version 3.0 onward replaces "i3-wm" with "wm" in Xresource keys.  The content on this page has been updated for Regolith 3.0+.  If you are using an earlier version, use "i3-wm" instead of "wm" in the key names below.  For example `wm.foo.bar` changes to `i3-wm.foo.bar` for Regolith 1.x and 2.x.
{{< /callout >}}

Die häufigste Änderung der Tastaturbelegung ist die Taste {{< keys "super" >}}. Regolith verwendet `Xresources` als
kanonische Quelle für Einstellungen, die von verschiedenen UI-Komponenten gelesen werden. Die Tabelle der `Xresources`
-Parameter, für die Benutzerkonfiguration [ist an anderer Stelle verfügbar]({{< ref "xresources" >}}). Um die
Tastaturbelegung {{< keys "super" >}} von der Taste "windows" auf "alt" zu ändern, fügen Sie die folgende Zeilen in die
Datei `~/.config/regolith3/Xresources` ein:

```yaml {filename="~/.config/regolith3/Xresources"}
wm.mod: Mod1
wm.alt: Mod4
```

**Anmerkung**: GNOME hat auch seinen eigenen Satz von Tastaturbelegungen. Wenn die Regolith-Sitzung zum ersten Mal
initialisiert wird, werden die kollidierenden GNOME-Tastaturbelegungen aus den Benutzereinstellungen entfernt.
GNOME-Tastaturbelegungen können in der Regolith-Einstellungs-App verwaltet werden.

![](/images/v-tour/regolith-gnome-keybindings.png)

## Systemverwaltung

Die App `regolith-control-center` ist das Werkzeug zum Konfigurieren von Gebietsschema, Datum, Anzeigen, Netzwerk und
verschiedenen anderen Einstellungen. Starten Sie sie über den App-Launcher mit {{< keys "super,space" >}}, geben
Sie `settings` ein und drücken Sie Enter, um die App zu starten. Die direkte Tastaturbelegung ist {{< keys "super,c" >}}. Die Anwendung kann auch von einem Terminal mit dem Befehl `regolith-control-center` gestartet werden.

![](/images/v-tour/regolith-settings-about.png)

## Weitere Lektüre

Wenn Sie tiefer einsteigen wollen, finden Sie die verfügbaren [Howtos]({{< ref "howtos" >}}) oder lesen Sie
die [Xresources reference]({{< ref "docs/Reference/xresources.md" >}}). Werden Sie ein i3-Poweruser, indem Sie
das [Benutzerhandbuch](https://i3wm.org/docs/userguide.html) lesen.
