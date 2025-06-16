---
title: 'i3設定のカスタマイズ'
linkTitle: 'i3設定のカスタマイズ'
description: >
    i3の外観や動作を変更する方法。
---

既定の設定は`/usr/share/regolith/i3/config.d`に保存されています。これらはアルファベット順に読み込まれます。次に、`~/.config/regolith3/i3/config.d`のユーザー設定がアルファベット順に読み込みまれます。Regolithの既定の設定は`Xresources`値の設定や、ユーザー設定の追加や、`apt`による既定の設定の追加・削除でカスタマイズできるように構築されています。これにより、Regolithではi3のカスタマイズ方法に、設定をクリーンにするため、個別に、もしくは組み合わせて使用できる複数のアプローチがあります。

-   `Xresources`によるキーバインド・文字列・プログラム名・色といった、オーバーライド値の使用
-   `apt`によるRegolithの既定ファイルを追加や削除
-   独自の設定ファイルの追加
-   Regolithの設定の大規模の置き換え

## `Xresources`変数を設定する方法

次の設定例（`/usr/share/regolith/i3/config.d/80_compositor`）より、`$wm.program.compositor`を設定することで、既存のi3設定ファイルを編集したりオーバーライドすることなく、使用したいコンポジターを読み込むようにカスタマイズすることができます。

{{< callout type="warning" >}}
Regolithバージョン3.0ではXresourceキーが"i3-wm"から"wm"に置き換えられます。このページの内容はRegolith 3.0以降のものに更新されています。それ以前のバージョンを使用する場合は、キーの名前を"wm"から"i3-wm"として使用してください。例として、Regolith 1.xと2.xでは`wm.foo.bar`は`i3-wm.foo.bar`に変更されます。
{{< /callout >}}

```text {filename="/usr/share/regolith/i3/config.d/80_compositor"}
set_from_resource $wm.program.compositor wm.program.compositor /usr/share/regolith-compositor/init
exec_always --no-startup-id $wm.program.compositor
```

これを行うには、`~/.config/regolith3/Xresources`に以下の行を追加するだけでできます。

```yaml {filename="~/.config/regolith3/Xresources"}
# 自作のコンポジターを使用する
wm.program.compositor: /usr/local/bin/my-compositor
```

同様に、`/usr/share/regolith/i3/config.d`のファイルを基にした`Xresources`で他のオプションのオーバーライドをすることもできます。

```yaml {filename="~/.config/regolith3/Xresources"}
## 間隔の設定
wm.gaps.inner.size: 1

## 分割線の設定
wm.window.border.size: 3
wm.client.focused.color.child_border: #AAD3E9

## ワークスペースの名前
wm.workspace.01.name: 1:FOO
```

詳しくは、[`Xresources`に関するより詳細なディスカッション]({{< ref "docs/howtos/override-xres.md" >}})を参照してください。

## 既定の部分設定をアンインストールする方法

提供された既定の部分設定を全体的に置き換えたり削除したい場合、提供された部分を`apt`でパッケージとして削除することができます。

削除したい部分設定を提供するパッケージの名前を特定するには次のコマンドを実行します。

```bash
dpkg -S /usr/share/doc/regolith-i3-workspace-config
```

出力結果にパッケージの名前が含まれている場合：

```yaml
regolith-i3-workspace-config: /usr/share/doc/regolith-i3-workspace-config
```

このパッケージで提供されるファイルの完全なリストを確認し、何を削除しようとしているのかわかっている場合：

```bash
dpkg -L regolith-i3-workspace-config
```

部分設定を削除するには：

```bash
sudo apt remove regolith-i3-workspace-config
```

## 既定の部分設定をカスタマイズする方法

{{< callout type="warning" >}}
Regolithの既定部分のみをカスタマイズする場合、`~/.config/regolith3/i3/config`ファイルが存在**しない**ことを確認してください。これはRegolithの既定をオーバーライドするため、このハウツーの内容と互換性がありません。
{{< /callout >}}

初めに、ユーザーの部分設定のディレクトリが存在することを確認します。

```bash
mkdir -p ~/.config/regolith3/i3/config.d
```

次に、既定の設定を独自の設定ディレクトリにコピーします。

```bash
cp /usr/share/regolith/i3/config.d/40_workspace-config ~/.config/regolith3/i3/config.d/
```

希望に沿った変更を加えます。

```bash
vim ~/.config/regolith3/i3/config.d/40_workspace-config
```

既定に含まれているパッケージをクリーンアップします。

```bash
sudo apt remove regolith-i3-workspace-config
```

最後に、**i3を再起動するか、ログアウトしてから戻ります**。

## Xresourcesオプションの範囲を超えてRegolithの既定を変更する方法

上記の記事の内容以上の設定が必要な場合は、[設定ファイルのステージング]({{< ref "stage-configs" >}})をご覧ください。

## Regolithの既定を完全に置き換える方法

Regolithの既定を保持したくない（例えば、あなたが完全な個人設定を使いたい経験豊富な`i3`ユーザーである）場合、すべてのRegolithの部分設定をアンインストールしたり、独自のRoot設定を`~/.config/regolith3/i3/config`に作成することができます。
