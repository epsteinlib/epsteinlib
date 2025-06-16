---
title: "設定ファイル"
linkTitle: "設定ファイル"
weight: 2
description: >
  Regolith使用が使用する設定ファイルについての一般的な情報
---

## 設定ファイルのクイックリファレンス

| **コンポーネント**              | **既定の設定**                                           | **ユーザー設定**                       | **注意**                                                                                                                                                                                             |
| -------------------------- | ------------------------------------------------------------ | ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| i3 / i3-gaps               | `/etc/regolith/i3/config`, `/usr/share/regolith/i3/config.d` | `~/.config/regolith3/i3`              | ユーザー設定ディレクトリで`config.d`部分設定を作成してRegolith2の設定を _拡張_ するか、`config`を作成して _全体にオーバーライド_ します。[詳細]({{< ref "docs/howtos/customize-i3-configuration.md" >}}) |
| Xresources                 | `/usr/share/regolith-look/default`                           | `~/.config/regolith3/Xresources`      | `~/.Xresources`も読み込まれますが、他のデスクトップセッションでも必要になる可能性のあるプロパティを対象にしています。[詳細]({{< ref "docs/howtos/override-xres.md" >}})                                    |
| バーステータスインジケーター      | `/usr/share/i3xrocks/conf.d/`                                | `~/.config/regolith3/i3xrocks/conf.d` | 各ブロックには独自の設定ファイルがありません。ファイル名のアルファベット順に読み込まれた設定です。[詳細]({{< ref "docs/howtos/add-remove-blocklets.md" >}})                                                                |
| Picom（Compton）コンポジター | `/etc/regolith/picom/config`                                 | `~/.config/regolith3/picom/config`    | [詳細]({{< ref "docs/howtos/customize-compositor.md" >}})                                                                                                                                             |

{{< callout type="warning" >}}
`~/.config/regolith3/i3/config`にi3の設定ファイルを追加すると、Xresourcesやlookを含むすべてのRegolithの設定がオーバーライドされます。
{{< /callout >}}

## 履歴

### 部分設定

Regolith 2.0では、`i3`が[4.20](https://i3wm.org/downloads/RELEASE-NOTES-4.20.txt)になったことで部分設定のサポートが追加されます。これは[includes](https://i3wm.org/docs/userguide.html#include)によって行われます。このincludesのサポートにより、Regolith 2.0がパッケージマネージャーを使用して、既定の部分設定を追加もしくは削除することができます。パッケージのインストールや削除により、 i3設定を特定の設定にカスタマイズしながら、設定を変更する必要のない部分については、上流の変更を追跡できます。

`/etc/regolith/i3/config`にあるrootのRegolith設定ファイルは、部分設定を次のように読み込みます：

```text {filename="/etc/regolith/i3/config"}
# Regolithのi3のあらゆる部分設定を読み込む
include /usr/share/regolith/i3/config.d/*

# ユーザーのi3の部分設定を読み込む
include $HOME/.config/regolith2/i3/config.d/*
```
