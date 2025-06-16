---
title: "バーステータスインジケーターの変更"
linkTitle: "バーステータスインジケーターの変更"
description: >
  バーにあるステータスインジケーターを変更する。
---

## ステータスインジケーターの追加と削除

### ステータスインジケーターの検索

```bash
apt search i3xrocks-
```

### ステータスインジケーターのインストール

```bash
sudo apt install i3xrocks-memory
regolith-look refresh
```

## i3バーのカスタマイズ {#customize-bar}

バーにある各インジケーターはファイルで管理されます。バーのステータスインジケーターの位置に対応しています。例として、ファイル名が`10_`で始まるものは、 `20_`で始まる設定ファイルのインジケーターより前に配置されることを意味します。各ファイルの中には、`i3xrocks`がスクリプトを実行するために使用する情報があり、最終的にバーに表示されるデータを返します。いくつかのステータスインジケーターは設定パラメーターが存在し、ユーザーの好みに合わせることができます。例えば、バッテリーの状態を頻繁に変更するのがよい場合、調べる間隔を更新することができます。バーをカスタマイズするために最初にするべきことは、表示したいインジケーターを`/usr/share/i3xrocks/conf.d/`から`~/.config/regolith3/i3xrocks/conf.d`にコピーすることです。 `~/.config/regolith3/i3xrocks/conf.d`内の各ファイルは好きなように変更できます。バーのステータスインジケーターを変更するには、ファイルの名前を好きな並び順に変更するだけです。

変更を加えたあとは、セッションを再読み込みすることで、新しい設定ファイルに基づいてバーが更新されるはずです。

例えば、バーにおけるバッテリーとトラフィックのブロックの変更と通知を非表示にしたいとき、以下の操作を行います。

```console
$ ls /usr/share/i3xrocks/conf.d/
01_setup
30_net-traffic
80_battery
80_rofication
90_time
$ mkdir -p ~/.config/regolith3/i3xrocks/conf.d
$ cp /usr/share/i3xrocks/conf.d/01_setup ~/.config/regolith3/i3xrocks/conf.d/01_setup
$ cp /usr/share/i3xrocks/conf.d/80_battery ~/.config/regolith3/i3xrocks/conf.d/30_battery
$ cp /usr/share/i3xrocks/conf.d/30_net-traffic ~/.config/regolith3/i3xrocks/conf.d/80_net-traffic
$ regolith-look refresh
```

上記の操作では、以下の操作が実行されます。

-   Regolithが作成したディレクトリに存在する現在の`i3xrocks`モジュールが一覧表示されます
-   ユーザー設定が存在するディレクトリに`~/.config/regolith3/i3xrocks/conf.d`が作成されます
-   3つのファイルが、アプリケーションが作成したディレクトリからユーザーが作成したディレクトリにコピーされます。これにより、この記事の冒頭で説明したようにユーザーが`n_`の番号を変更すると、順番が変更されます。
-   `regolith-look refresh`コマンドで**バーステータスインジケーター**の再読み込みが行われます

**注意**：ユーザーディレクトリの`~/.config/regolith3/i3xrocks/conf.d`に複数のブロック設定が存在する場合、既定である`/usr/share/i3xrocks/conf.d/`は無視されます。

## さらに深く読み込む

[設定のリファレンスページ]({{< ref "/docs/Reference/configurations.md" >}})には、Regolith Linuxを使うときの設定ファイルについての詳細な項目があります。
