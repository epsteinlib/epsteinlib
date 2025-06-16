---
title: "Grundlegende Verwendung"
linkTitle: "Grundlegende Verwendung"
weight: 3
description: >
 Lernen Sie zunächst die wichtigsten Teile von Regolith kennen.
---

Während wir in anderen Bereichen der Dokumentation alles detaillierter behandeln werden, deckt diese Seite die grundlegende Verwendung von Regolith ab.

## Starten einer Anwendung

Nach der ersten Anmeldung enthält der Bildschirm eine Leiste am unteren Rand und einige Tastaturkürzel in einem Dialog. Dieser Dialog soll dem Benutzer helfen, Tastenkombinationen für verschiedene Aktionen in der Desktop-Umgebung zu finden und zu speichern. Er kann mit der Taste `esc` verlassen und mit {{< keys "super,shift,?" >}}. Um eine Anwendung zu
starten, drücken Sie {{< keys "super,space" >}} und ein ähnlicher Dialog wird mit den auf dem System installierten Anwendungen geladen. Wenn Sie beginnen, den Namen der Anwendung, die Sie starten möchten, einzugeben, wird die Liste auf die passenden Anwendungen eingeschränkt. Sie können auch mit den Pfeiltasten durch die Liste blättern.

![](/images/v-tour/regolith-ilia-apps.png)

## Zu einem Fenster navigieren

Um eine laufende Anwendung zu fokussieren, kann die {{< keys "super,ctrl,space" >}} Tastenkombination verwendet werden, um einen Dialog mit allen laufenden Fenstern anzuzeigen. Wählen Sie einfach das gewünschte Fenster aus, und das System wechselt in den Arbeitsbereich des Fensters und aktiviert dessen Fokus.

![](/images/v-tour/regolith-ilia-windows.png)

## Einen neuen Arbeitsbereich erstellen

Arbeitsbereiche (in manchen Desktop-Umgebungen auch als "virtuelle Desktops" bezeichnet) sind eine der wichtigsten Möglichkeiten, um Fenster zu verwalten. Jeder Arbeitsbereich bietet einen Bildschirm für einen unabhängigen Satz von Fenstern. Um einen Arbeitsbereich zu erstellen oder zu einem bestehenden Arbeitsbereich zu navigieren, geben Sie seinen
Index mit den Tasten {{< keys "super,0-9" >}} an.

## Suchen und Starten einer Datei

Die Suche nach einer bestimmten Datei auf Ihrem lokalen Computer kann auf verschiedene Weise erfolgen. Regolith bietet einen Dateisuchdialog über {{< keys "super,alt,space" >}} zum schnellen Suchen und optionalen Laden einer Datei mit seiner Standardanwendung. Durch Drücken der Eingabetaste bei einer Auswahl wird die Datei mit dem Befehl `xdg-open` aufgerufen.

![](/images/v-tour/regolith-ilia-files.png)

## Terminal

Das Terminal ist eine spezielle Anwendung, die ihre eigene Tastenkombination hat, {{< keys "super,enter" >}}.

## Browser

Der bevorzugte Browser ist ebenfalls eine spezielle Anwendung und kann mit {{< keys "super,shift,enter" >}} gestartet werden.

## Beenden einer Anwendung

Mit {{< keys "super,shift,q" >}} wird die fokussierte Anwendung beendet. Bei fehlerhaften Anwendungen wird {{<
keys "super,alt,q" >}} das Beenden der fokussierten Anwendung erzwingen, aber beachten Sie, dass ungespeicherte Daten
verloren gehen können.

## Arbeitsbereiche wechseln

Eine primäre Möglichkeit, die Belegung des Bildschirms durch Anwendungen zu verwalten, ist die Verwendung von Arbeitsbereichen (Workspaces). Ein Arbeitsbereich ist im Wesentlichen eine Gruppe von null oder mehr Anwendungen, die alle Fenster auf demselben virtuellen Bildschirm darstellen. Benutzer können zwischen Arbeitsbereichen wechseln, um mehr
Platz für Anwendungen zu gewinnen. Fenster können zwischen Arbeitsbereichen verschoben werden. Arbeitsbereiche können an eine physische Anzeige gebunden sein oder nach Bedarf auf andere Anzeigen verschoben werden. Um den Arbeitsbereich zu wechseln, halten Sie die Taste {{< keys "super" >}} gedrückt und drücken Sie eine der Zifferntasten. Dadurch wird zu diesem Arbeitsbereich navigiert. Wenn in einem bestimmten Arbeitsbereich ein Fenster vorhanden ist, wird es in der Leiste als aktiv angezeigt. So kann der Benutzer anhand der Leiste erkennen, in welchen Arbeitsbereichen sich aktive Anwendungen befinden.

Arbeitsbereiche können auch mit {{< keys "super,tab" >}} und {{< keys "super,shift,tab" >}} vor- und zurückgeblättert werden.

## Das System konfigurieren

Die Einstellungs-App (auch bekannt als `regolith-control-center`) ist die primäre Schnittstelle für allgemeine Konfigurationsaufgaben, wie die Konfiguration von WiFi und Bluetooth, die Konfiguration von App-Zuordnungen und die Einrichtung von Benutzern. Die App kann über die {{< keys "super,c" >}} Tastenkombination oder vom Launcher aus über die App `Regolith Settings` gestartet werden.

![](/images/v-tour/regolith-gnome-settings.png)

## Tastenkombinationen anzeigen

Alle vom Regolith-Desktop verwendeten Tastenkombinationen können mit `ilia` eingesehen, sortiert und ausgeführt werden. Dieser Dialog wird mit {{< keys "super,shift,?" >}} aktiviert.

![](/images/v-tour/regolith-ilia-keybindings.png)

## Verwalten von Benachrichtigungen

Das Desktop-Benachrichtigungssystem von Regolith wurde entwickelt, um Pop-ups zu verhindern, die den Benutzer von seiner Arbeit ablenken könnten. Benachrichtigungen können über {{< keys "super,n" >}} eingesehen und verwaltet werden. Spezielle Tastenkombinationen, wie z.B. das Löschen einer oder mehrerer Benachrichtigungen, finden Sie auf der Registerkarte "Hilfe" des Dialogs. Für Benutzer, die eine bessere Sichtbarkeit der aktuellen Benachrichtigungen wünschen, kann der Statusindikator `i3xrocks-rofication` über den Systempaketmanager installiert werden. Damit wird die Anzahl der kritischen und normalen ungelesenen Benachrichtigungen in der Leiste angezeigt.
