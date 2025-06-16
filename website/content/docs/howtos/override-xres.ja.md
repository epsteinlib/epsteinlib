---
title: "Xresourcesのオーバーライド"
linkTitle: "Xresourcesのオーバーライド"
description: >
  ユーザーがコピーしたRegolithの設定ファイルをステージする方法について学ぶ
---

Regolithは、[Xresourcesシステム](https://en.wikipedia.org/wiki/X_resources)に依存し、統合されたインターフェース設定を提供しています。`Xresources`の値を変更することで、Regolithは、ユーザーインターフェイスの更新、カスタム動作の指定、バーの一部のblockletの特定の形式の定義などの方法でカスタマイズできます。

{{< callout type="info" >}}
Regolithでは、`Xresources`の値は、独自の`Xresources`ファイルを用意する方法か、既存の`Xresources`の値をオーバーライドする方法の2つの方法で変更できます。後者のアプローチは、Regolithが新しいバージョンにアップグレードする際にカスタマイズを妨げる可能性を下げるため、強く推奨されます。`Xresources`ツリーを完全にコピーする場合、将来のバージョンで破損した変更を手動で統合する必要があります。
{{< /callout >}}

## 初期化

Regolithセッションが開始されるとき、以下の順番でXresourcesが読み込まれます：

| ファイル                            | オプション | 注意                                                                             |
| ------------------------------- | -------- | --------------------------------------------------------------------------------- |
| `~/.Xresources`                 | Y        | Regolith以外の設定を対象にしています                                                |
| `/usr/share/regolith-look/**/`     | N        | Regolithの`Xresources`は外観と関連付けられています |
| `~/.config/regolith3/Xresources` | Y        | 特定のオーバーライドを`Xresources`の既定に適用します                               |

冗長設定の指定を必要とせず、長期にわたって管理しやすいため、カスタマイズには`~/.config/regolith3/Xresources`を使用することを推奨します。

## 変更可能な値の決定

`xrdb`ツールは既存の`Xresources`の値をリスト表示するときに使用することができます。Regolith 2.1リリースには、[いくつか既存のオーバーライドが存在します]({{< ref "/docs/Reference/xresources.md" >}})：

```console
$ xrdb -query
[...]
gnome.terminal.scrollbar:	never
gnome.terminal.use-transparent-background:	true
gtk.document_font_name:	Sans 12
gtk.font_name:	Bitstream Vera Sans 12
gtk.icon_theme_name:	Papirus-Dark
gtk.monospace_font_name:	BitstreamVeraSansMono Nerd Font 13
gtk.theme_name:	Ayu-Mirage-Dark
wm.bar.background.color:	#1F2430
[...]
```

## 参考例

以下のコマンドはファイルの最後にテキストを付け加えるため、コマンドを複数回実行すると、重複した入力が発生することに注意してください。

### 例 - UIを高DPI画面に更新する

`~/.config/regolith3/Xresources`のオーバーライドファイルを使用する場合、変更したい値を指定するだけです。`xrdb`ツールを使用して、現在の値を変更することもできます。

1. `~/.config/regolith3/Xresources`ファイルに次の値を作成もしくは追加する：

```yaml {filename="~/.config/regolith3/Xresources"}
Xft.dpi: 192
```

2. `Xresources`設定を再読込する:

```bash
regolith-look refresh
```

3. 新規のターミナルを開き、変更が有効になっていることを確認する。

`192`は参考例の値であり、必要に応じて変更してください。

### 例 - i3バーの位置を変更する

```console
$ xrdb -query | grep position
wm.bar.position:	bottom
$ echo "wm.bar.position:	top" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### 例 - GTKテーマを変更する

```console
$ xrdb -query | grep gtk
gtk.theme_name:	Ayu-Mirage-Dark
$ echo "gnome.gtk.theme:	Adwaita" >> ~/.config/regolith3/Xresources
$ regolith-look refresh
```

### 例 - システムトレイを無効にする

```bash
echo "wm.bar.trayoutput:	none" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### 例 - AltキーをWinキーの代わりにSuperキーとして使用する

```bash
echo "wm.mod: Mod1" >> ~/.config/regolith3/Xresources
echo "wm.alt: Mod4" >> ~/.config/regolith3/Xresources
```

次に、変更を有効にするため、i3を再読み込みします。i3が別のプロセスとして起動するものなど、一部の設定では、セッションに再ログインする必要があります。

## さらに深く読み込む

Regolithの設定ファイルについて詳細については[設定のリファレンスページ]({{< ref "/docs/Reference/configurations.md" >}})を参照してください。
