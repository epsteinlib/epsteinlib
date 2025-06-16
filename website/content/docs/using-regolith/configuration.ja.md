---
title: '設定'
linkTitle: '設定'
weight: 4
description: >
    Regolithの外観と動作を変更する。
next: /docs/howtos
---

## 新規ユーザーに推奨されるパッケージ {#recommended-packages-for-new-users}

### 外観

Regolith公式の外観をすべてインストールする:

```bash
sudo apt install regolith-look-*
```

もしくは、`apt`でインストールされて利用可能な外観のリストを表示する:

```bash
apt list | grep regolith-look-
```

注意：これはフォントやテーマ、壁紙も含みます。ディスクの使用率の最小化を望むユーザーは、使用していない外観をアンインストールすることを推奨します。

外観の選択方法の概要については、以下の項目を参照してください。

### ステータスインジケーター

次の行では、新規ユーザーに推奨されるステータスインジケータパッケージをインストールします。

```bash
sudo apt install i3xrocks-focused-window-name i3xrocks-rofication i3xrocks-info i3xrocks-app-launcher i3xrocks-memory
```

これらのインジケーターを後から削除したいときは、システムからアンインストールします。例：

```bash
sudo apt remove i3xrocks-info
```

#### ラップトップのバッテリー状態

ラップトップの現在のバッテリー状態を表示するには：

```bash
sudo apt install i3xrocks-battery
```

### コンポジター

コンポジターはウィンドウ環境に仮想エフェクトを追加します。既定では、Regolith 3には`picom`コンポジターがインストールされています。特定のハードウェアやドライバの制限によってグラフィックの問題が発生した場合、変わりに次のコンポジターをインストールする必要があります。

| パッケージ                      | 注意          |
|------------------------------|-------------------|
| `regolith-compositor-picom-glx` | 既定 |
| `regolith-compositor-compton-glx` | 古いコンポジターであり、いくつかのハードウェアでよりよく動作する場合があります |
| `regolith-compositor-xcompmgr` | 古いコンポジターであり、いくつかのハードウェアでよりよく動作する場合があります |
| `regolith-compositor-none` | パススルーパッケージであり、コンポジターは実行されません |

詳細な概要は、[このハウツー]({{< ref "docs/howtos/customize-compositor.md" >}})を参照してください。

## 壁紙

多くのRegolithの外観は、既定のデスクトップカラーと画像を提供します。壁紙はXressourcesキー・`regolith.wallpaper.file`に壁紙画像のパスを指定することで設定することができます。  

例として、`/usr/share/backgrounds/hardy_wallpaper_uhd.png`ファイルが存在すると仮定したとき： 

```bash
echo "regolith.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

上記のコマンドでは、前述のプロパティを使用するときに他のパスから`/usr/share/backgrounds/hardy_wallpaper_uhd.png`に変更されます。

壁紙画像の表示形式（拡大、大きさ、など…）を変更するには、`regolith.wallpaper.options`を設定します。オプション：`none`, `wallpaper`, `centered`, `scaled`, `stretched`, `zoom`, `spanned`

例：

```bash
echo "regolith.wallpaper.options: zoom" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

Xresourceキー・`regolith.wallpaper.color.primary`を使用して、画像ではなく色を指定するには：

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

2つ目の色とグラデーションを指定するには：

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: blue" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.secondary: green" >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.shading.type: vertical" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### ロック画面の壁紙

ロック画面壁紙はデスクトップ壁紙と同じ方法で管理することができます。キーはデスクトップ壁紙と似たようなものですが、`regolith.lockscreen.wallpaper...`のように、`lockscreen`を`wallpaper`の前に追加します。例えば、ロック画面の画像を指定する場合：

```bash
echo "regolith.lockscreen.wallpaper.file: /usr/share/backgrounds/hardy_wallpaper_uhd.png" >> ~/.config/regolith3/Xresources
regolith-look refresh
```

### 壁紙の処理を無効化

Regolithの外部で壁紙を管理したい場合、壁紙画像と色の値を空にします。

```bash
echo "regolith.wallpaper.file: " >> ~/.config/regolith3/Xresources
echo "regolith.wallpaper.color.primary: " >> ~/.config/regolith3/Xresources
```

## ステータスバーのインジケーター

CPUの読み込み、日付、時間、通知、天気や他のシステム情報などのステータスインジケーターは、パッケージのインストールで追加や削除をすることが可能です。例えば、ラップトップバッテリーのステータスインジケーターを表示するには、`sudo apt install i3xrocks-battery`を実行して、セッションを`regolith-look refresh`で再読み込みするだけです。利用可能なインジケーターを探すには、`apt search ^i3xrocks-`を実行するか、[Synaptic](https://help.ubuntu.com/community/SynapticHowto)のようなGUIのパッケージマネージャーで`i3xrocks-`と検索して見つけることができます。「ブロックレット(blocklet)」と呼ばれる、ステータスインジケーターを設定するのに[利用可能なドキュメント]({{< ref "/docs/howtos/add-remove-blocklets.md" >}})が存在します。以下は、Regolith 2.2時点でのステータスインジケーターのリストです。

| パッケージ                      | ステータスの機能 | 
|--------------------------------|-------------------|
| i3xrocks-app-launcher          | アプリ起動ダイアログのショートカット |                               
| i3xrocks-battery               | バッテリーの状態 |               
| i3xrocks-bluetooth             | Bluetoothの状態 |                             
| i3xrocks-cpu-usage             | CPU使用率 |                             
| i3xrocks-disk-capacity         | ディスク使用率 |                                 
| i3xrocks-focused-window-name   | フォーカスされたウィンドウのタイトル |                                       
| i3xrocks-info                  | キーバインドビューワーのショートカット |            
| i3xrocks-key-indicator         | キーボード切り替えの状態 |                                 
| i3xrocks-keyboard-layout       | キーボードレイアウト |                                   
| i3xrocks-media-player          | メディアプレーヤーの状態 |                               
| i3xrocks-memory                | メモリ使用率 |              
| i3xrocks-microphone            | マイクの状態 |                             
| i3xrocks-net-traffic           | ネットワークの状態 |                               
| i3xrocks-next-workspace        | 新規ワークスペース作成のショートカット |                                 
| i3xrocks-nm-vpn                | VPNの状態 |                         
| i3xrocks-openvpn               | OpenVPNの状態 |               
| i3xrocks-rofication            | 通知の状態 |                             
| i3xrocks-tailscale             | Tailscale status |                             
| i3xrocks-temp                  | システムの温度 |                       
| i3xrocks-time                  | システム時間 |                       
| i3xrocks-updates               | システムパッケージの状態 |                           
| i3xrocks-volume                | 音量の状態 |                         
| i3xrocks-weather               | 地域の天気 |                           
| i3xrocks-wifi                  | 接続中のWiFiアクセスポイント名 |                       

## 外観

色、壁紙、ウィンドウとバーのレイアウトや他の視覚的要素はRegolithにバンドルされ、"外観(look)"と呼ばれています。外観は、デスクトップの全体の外観を変更する簡単な方法を提供します。ステータスバーのインジケーターのような外観はパッケージ化され、他のソフトウェアのようにインストールしたりアンインストールしたりすることができます。外観パッケージは命名フォーマットに従い`regolith-look-<name>`となっています。`apt`やGUIパッケージマネージャーから利用可能な外観を検索することができます。{{< keys "super,alt,l" >}}キーバインドで、ダイアログからインストールされた外観を選択します。

![](/images/v-tour/regolith-ilia-look-selector.png)

また、ターミナルのコマンド・`regolith-look`で外観を変更したり、選択した外観でアクティブなセッションを再読み込みすることができます。ここでは例として、`gruvbox`の外観に切り替えます。

```console
$ apt search ^regolith-look-
[...]
regolith-look-gruvbox/unknown,now 0.4.6-1regolith amd64
[...]
$ sudo apt install regolith-look-gruvbox
regolith-look set gruvbox
regolith-look refresh
```

利用可能な外観をすべてインストールする簡単な方法：

```bash
sudo apt install regolith-look-*
```

## i3の機能

[Regolith 2.0の開始時から]({{< ref "docs/reference/configurations.md#history" >}})、いくつかのi3の外観設定はパッケージマネージャーによって管理されています。Regolithでは、i3の設定がすべて「部分設定(config partial)」として提供されます。これらのパッケージはi3設定の一部として、`/usr/share/regolith/i3/config.d`以下にインストールされます。パッケージのインストールと削除をすることで、i3設定を特定の設定に合わせてカスタマイズすることができます。既定では、 `regolith-desktop`パッケージがインストールされたときに設定要素が依存関係としてインストールされています。

| パッケージ                      | 機能          | ファイル |
|------------------------------|-------------------|----- |
| regolith-i3-base-launchers   | ブラウザとターミナルの起動 | `/usr/share/regolith/i3/config.d/15_base_launchers` |
| regolith-i3-default-style    | ウィンドウの振る舞い | `/usr/share/regolith/i3/config.d/40_default-style` |
| regolith-i3-navigation       | 操作のキーバインド | `/usr/share/regolith/i3/config.d/30_navigation` |
| regolith-i3-networkmanager   | ネットワークとWiFiの機能 | `/usr/share/regolith/i3/config.d/88_network-manager` |
| regolith-i3-next-workspace   | ワークスペースのユーティリティー | `/usr/share/regolith/i3/config.d/40_next-workspace` |
| regolith-i3-resize           | リサイズのキーバインド | `/usr/share/regolith/i3/config.d/50_resize-mode` |
| regolith-i3-session          | セッションのキーバインド | `/usr/share/regolith/i3/config.d/55_session_keybindings` |
| regolith-i3-unclutter        | マウスカーソルを使用していないときに非表示化 | `/usr/share/regolith/i3/config.d/86_unclutter` |
| regolith-i3-user-programs    | オプションでXresourcesで指定されたユーザープログラムの起動 | `/usr/share/regolith/i3/config.d/90_user-programs` |
| regolith-i3-workspace-config | ワークスペースのキーバインド | `/usr/share/regolith/i3/config.d/40_workspace-config` |

ソフトの依存関係は、それに依存するパッケージを削除することなく削除できます。これはリストにあるいくつかのパッケージを削除できることを示しています。ユーザーは、デフォルトの設定部分を独自のバージョンに置き換えることで、より安定した方法で設定を微調整できます。これはカスタマイズしたい部分を`~/.config/regolith3/i3/config.d/`にコピーして、`apt`で`/usr/share/regolith/i3/config.d`のオリジナルのファイルを削除することでできます。この方法により、変更されていない設定部分は、ユーザー固有の設定に影響を与えることなく更新やバグ修正を取得し続けることができます。

参考例は[i3設定のカスタマイズのハウツー]({{< ref "docs/howtos/customize-i3-configuration.md" >}})を参照してください。


### すべてのi3設定パッケージ

ここに記載されているすべてのi3設定パッケージはRegolith 2.1で利用できます。

| パッケージ                      | 機能          | ファイル |
|------------------------------|-------------------|----- |
| regolith-i3-base-launchers   | ターミナルとブラウザを起動 | `/usr/share/regolith/i3/config.d/15_base_launchers` |
| regolith-i3-compositor       | コンポジターの統合（`regolith-compositor-<variant>`へ転送） | `/usr/share/regolith/i3/config.d/80_compositor` |
| regolith-i3-default-style    | ウィンドウの振る舞い | `/usr/share/regolith/i3/config.d/40_default-style` |
| regolith-i3-ftue             | 最初のユーザーエクスペリエンス | `/usr/share/regolith/i3/config.d/84_ftue` |
| regolith-i3-gaps-partial     | `i3-gaps`設定 （インストールした場合i3がi3-gapsに置き換えられます） | `/usr/share/regolith/i3/config.d/35_i3-gaps` |
| regolith-i3-gnome            | Regolith設定のキーバインド | `/usr/share/regolith/i3/config.d/60_config_keybindings` |
| regolith-i3-i3xrocks         | i3バーのステータスインジケーターの設定 | `/usr/share/regolith/i3/config.d/70_bar` |
| regolith-i3-ilia             | `ilia`の統合 | `/usr/share/regolith/i3/config.d/20_ilia` |
| regolith-i3-navigation       | 操作のキーバインド | `/usr/share/regolith/i3/config.d/30_navigation` |
| regolith-i3-networkmanager   | ネットワークとWiFiの機能 | `/usr/share/regolith/i3/config.d/88_network-manager` |
| regolith-i3-next-workspace   | ワークスペースのユーティリティー | `/usr/share/regolith/i3/config.d/40_next-workspace` |
| regolith-i3-remontoire       | 古いキーバインドビューとの統合 | `/usr/share/regolith/i3/config.d/45_remontoire` |
| regolith-i3-resize           | リサイズのキーバインド | `/usr/share/regolith/i3/config.d/50_resize-mode` |
| regolith-i3-rofi             | 古いデスクトップランチャーとの統合 | `/usr/share/regolith/i3/config.d/20_rofi` |
| regolith-i3-rofication       | 古い通知ビューワーとの統合 | `/usr/share/regolith/i3/config.d/82_rofication` |
| regolith-i3-rofication-ilia  | 通知ビューワーとの統合 | `/usr/share/regolith/i3/config.d/82_rofication-ilia` |
| regolith-i3-session          | セッションのキーバインド | `/usr/share/regolith/i3/config.d/55_session_keybindings` |
| regolith-i3-snapshot         | i3ワークスペースの固定 | `/usr/share/regolith/i3/config.d/75_i3-snapshot` |
| regolith-i3-swap-focus       | 最後にフォーカスしたウィンドウへの切り替え | `/usr/share/regolith/i3/config.d/40_i3-swap-focus` |
| regolith-i3-unclutter        | マウスを使用していないときに非表示化 | `/usr/share/regolith/i3/config.d/86_unclutter` |
| regolith-i3-user-programs    | オプションでXresourcesで指定されたユーザープログラムの起動 | `/usr/share/regolith/i3/config.d/90_user-programs` |
| regolith-i3-workspace-config | ワークスペースのキーバインド | `/usr/share/regolith/i3/config.d/40_workspace-config` |

## キーバインド

{{< callout type="warning" >}}
Regolithバージョン3.0ではXresourceキーが"i3-wm"から"wm"に置き換えられます。このページの内容はRegolith 3.0以降のものに更新されています。それ以前のバージョンを使用する場合は、キーの名前を"wm"から"i3-wm"として使用してください。例として、Regolith 1.xと2.xでは`wm.foo.bar`は`i3-wm.foo.bar`に変更されます。
{{< /callout >}}

もっとも重要なキーバインドの変更は{{< keys "super" >}}キーです。Regolithは、さまざまなUIコンポーネントで読み取られる設定の真の正規ソースとして`Xresources`を利用できます。ユーザー設定で利用可能な`Xresources`キーのリストは、[ここから参照することができます]({{< ref "xresources" >}})。既定の{{< keys "super" >}}バインドを"windows"キーから"alt"に変更するには`~/.config/regolith3/Xresources`に以下の行を追加します。

```yaml {filename="~/.config/regolith3/Xresources"}
wm.mod: Mod1
wm.alt: Mod4
```

**注意**：GNOMEには独自のキーバインドも存在します。Regolithセッションが最初に初期化されると、 競合するGNOMEキーバインドがユーザー設定から削除されます。GNOMEのキーバインドは、Regolithの設定アプリから管理することができます。

![](/images/v-tour/regolith-gnome-keybindings.png)

## システム管理

`regolith-control-center`アプリは、言語、日付、ディスプレイ、ネットワークや他のさまざまな種類の設定ツールです。 {{< keys "super,space" >}}でアプリランチャーを開いて`settings`と入力し、起動したいアプリをヒットさせることで起動します。{{< keys "super,c" >}}は直接実行できるキーバインドです。`regolith-control-center`とターミナルに入力することでもアプリケーションを起動することができます。

![](/images/v-tour/regolith-settings-about.png)

## さらに深く読み込む

もっと深く知りたい人は、[このハウツー]({{< ref "howtos" >}})を利用したり、[Xresourcesリファレンス]({{< ref "docs/Reference/xresources.md" >}})を読みましょう。 [ユーザーガイド（英語）](https://i3wm.org/docs/userguide.html)を読んでi3のパワーユーザーになりましょう。

もしXresourcesが提供する以上の高度な設定が必要な場合は、[設定ファイルのステージング]({{< ref "stage-configs" >}})で見つけたi3設定のハウツーを利用してください。
