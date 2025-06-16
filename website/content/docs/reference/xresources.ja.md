---
title: "Xresources"
linkTitle: "Xresources"
weight: 7
description: >
  Regolithにおけるすべての`Xresources`キーの定義
---

{{< callout type="warning" >}}
Regolithバージョン3.0ではXresourceキーが"i3-wm"から"wm"に置き換えられます。このページの内容はRegolith 3.0以降のものに更新されています。それ以前のバージョンを使用する場合は、キーの名前を"wm"から"i3-wm"として使用してください。例として、Regolith 1.xと2.xでは`wm.foo.bar`は`i3-wm.foo.bar`に変更されます。
{{< /callout >}}

## Main `Xresources` keys

| `Xresources`キー                                   | 概要                                                                                                         |
| -------------------------------------------------- | ------------------------------------------------------------------- |
| gtk.theme_name                                     | GTKテーマ |
| gtk.icon_theme_name                                | GTKアイコンテーマ  |
| gnome.terminal.audible-bell                        | gnome-terminalのベルを有効にする（`true / false`）                       |
| gtk.document_font_name                             | GTKドキュメントのフォント名 |
| gtk.font_name                                      | GTKアプリとデスクトップにおける一般的な用途のフォント |
| gtk.monospace_font_name                            | i3barやターミナルなどのモノスペースエリアのフォント |
| gnome.terminal.background-transparency-percent     | gnome-terminalの背景の透過度（`0 - 100`）。`gnome.terminal.use-transparent-background`も参照。 |
| gnome.terminal.font                                | gnome-terminalのフォント定義                                                                                  |
| gnome.terminal.scrollbar                           | スクロールバーの有効・無効フラグ（`never / always`）                                                            |
| gnome.terminal.use-transparent-background          | gnome-terminalの背景の透過を有効にする（`true / false`）。`gnome.terminal.background-transparency-percent`も参照。 |
| gnome.wallpaper                                    | GNOME壁紙（設定をオーバーライド可能)                                                                     |
| gnome.wm.theme                                     | GNOMEウィンドウマネージャーのテーマ（使用されていません）                                                                                 |
| wm.alt                                          | Altキーのキーマップ。[i3ユーザーガイド](https://i3wm.org/docs/userguide.html#_using_i3)の概要も参照。       |
| wm.bar.background.color                         | バーの背景色                                                                                                |
| wm.bar.position                                 | バーの表示位置                                                                                                 |
| wm.bar.separator.color                          | ブロック間の分割線の色                                                                                   |
| wm.bar.separatorchar                            | ブロック間に表示する字体                                                                                 |
| wm.bar.statusline.color                         |                                                                                                                     |
| wm.bar.trayoutput                               | i3barの[アプリトレイ](https://i3wm.org/docs/userguide.html#_tray_output)有効・無効                        |
| wm.bar.workspace.active.background.color        |                                                                                                                     |
| wm.bar.workspace.active.border.color            |                                                                                                                     |
| wm.bar.workspace.active.text.color              |                                                                                                                     |
| wm.bar.workspace.focused.background.color       |                                                                                                                     |
| wm.bar.workspace.focused.border.color           |                                                                                                                     |
| wm.bar.workspace.focused.text.color             |                                                                                                                     |
| wm.bar.workspace.inactive.background.color      |                                                                                                                     |
| wm.bar.workspace.inactive.border.color          |                                                                                                                     |
| wm.bar.workspace.inactive.text.color            |                                                                                                                     |
| wm.bar.workspace.urgent.background.color        |                                                                                                                     |
| wm.bar.workspace.urgent.border.color            |                                                                                                                     |
| wm.bar.workspace.urgent.text.color              |                                                                                                                     |
| wm.client.focused.color.background              |                                                                                                                     |
| wm.client.focused.color.border                  |                                                                                                                     |
| wm.client.focused.color.child_border            |                                                                                                                     |
| wm.client.focused.color.indicator               |                                                                                                                     |
| wm.client.focused.color.text                    |                                                                                                                     |
| wm.client.focused_inactive.color.background     |                                                                                                                     |
| wm.client.focused_inactive.color.border         |                                                                                                                     |
| wm.client.focused_inactive.color.child_border   |                                                                                                                     |
| wm.client.focused_inactive.color.indicator      |                                                                                                                     |
| wm.client.focused_inactive.color.text           |                                                                                                                     |
| wm.client.unfocused.color.background            |                                                                                                                     |
| wm.client.unfocused.color.border                |                                                                                                                     |
| wm.client.unfocused.color.child_border          |                                                                                                                     |
| wm.client.unfocused.color.indicator             |                                                                                                                     |
| wm.client.unfocused.color.text                  |                                                                                                                     |
| wm.client.urgent.color.background               |                                                                                                                     |
| wm.client.urgent.color.border                   |                                                                                                                     |
| wm.client.urgent.color.child_border             |                                                                                                                     |
| wm.client.urgent.color.indicator                |                                                                                                                     |
| wm.client.urgent.color.text                     |                                                                                                                     |
| wm.floatingwindow.border.size                   | 浮遊しているウィンドウの分割線のサイズ（ピクセル）                                                                   |
| wm.font                                         | ウィンドウマネージャーのフォント                                                                                                 |
| wm.gaps.inner.size                              | 既定のウィンドウ間の隙間（ピクセル）                                                                          |
| wm.gaps.outer.size                              | 既定のスクリーンの端との隙間（ピクセル）                                                                           |
| wm.mod                                          | Superキーのキーマップ。[i3ユーザーガイド](https://i3wm.org/docs/userguide.html#_using_i3)の概要も参照。     |
| wm.window.border.size                           | 既定の分割線サイズ（ピクセル）                                                                                       |
| wm.workspace.01.name .. wm.workspace.19.name | ワークスペースのラベル                                                                                                    |
| wm.workspace.01.key .. wm.workspace.10.key   | ワークスペースのショートカット                                                                                                 |
| i3xrocks.action.volume.left                        | インジケーター上で左クリックをしたときのコマンド（既定：なし）                                           |
| i3xrocks.action.volume.middle                      | インジケーター上で中クリックをしたときのコマンド（既定：なし）                                          |
| i3xrocks.action.volume.right                       | インジケーター上で右クリックをしたときのコマンド（既定：音量のミュート）                            |
| i3xrocks.action.volume.scrollup                    | 音量インジケーターで上方向にスクロールしたときのコマンド（既定：音量を上げる）       |
| i3xrocks.action.volume.scrolldn                    | 音量インジケーターで下方向にスクロールしたときのコマンド（既定：音量を下げる）     |
| i3xrocks.error.color                               | エラーモードにおけるバーモジュールの色                                                                                  |
| i3xrocks.label.battery0                            | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.battery100                          | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.battery20                           | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.battery50                           | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.battery80                           | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.batteryac                           | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.batterycharging                     | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.batteryunknown                      | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.capslock                            | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.color                               | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.cpu                                 | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.disk                                | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.dn                                  | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.help                                | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.mediapaused                         | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.mediaplaying                        | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.memory                              | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.mute                                | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.notify.error                        | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.notify.none                         | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.notify.some                         | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.numlock                             | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.off                                 | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.on                                  | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.sound                               | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.thermometer                         | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.time                                | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.up                                  | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.wifi                                | バーアイコンの絵文字                                                                                                      |
| i3xrocks.label.wired                               | バーアイコンの絵文字                                                                                                      |
| i3xrocks.nominal                                   | 既定のバーアイコンの色                                                                                              |
| i3xrocks.value.color                               | バーモジュールにおける数値の色                                                                              |
| i3xrocks.value.font                                | バーモジュールにおける数値のフォント                                                                              |
| i3xrocks.warning                                   | バーモジュールの警告色                                                                                        |

## オプションの`Xresources`キー

次の`Xresources`キーは既定では未定義ですが、ユーザーは必要に応じて動作を変更するために使用できます。

| `Xresources`キー                  | 概要                                                                                                                                                                                         |
| --------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| wm.gaps.focus_follows_mouse    | 有効であればマウスの下のウィンドウをフォーカスする                                                                                                                                                     |
| wm.program.launcher.app        | アプリランチャーキーバインドを実行するi3コマンド                                                                                                                                                |
| wm.program.launcher.cmd        | コマンドランチャーキーバインドを実行するi3コマンド                                                                                                                                            |
| wm.program.launcher.window     | ウィンドウ操作キーバインドを実行するi3コマンド                                                                                                                                            |
| wm.program.help                | ヘルプキーバインドを実行するコマンド                                                                                                                                                          |
| wm.program.file_search         | ファイル検索キーバインドを実行するコマンド                                                                                                                                                   |
| wm.program.refresh_ui          | セッションの再読み込みキーバインドを実行するコマンド                                                                                                                                               |
| wm.program.logout              | ログアウトキーバインドを実行するコマンド                                                                                                                                                        |
| wm.program.reboot              | 再起動キーバインドを実行するコマンド                                                                                                                                                        |
| wm.program.shutdown            | シャットダウンキーバインドを実行するコマンド                                                                                                                                                      |
| wm.program.lock                | 画面ロックキーバインドを実行するコマンド                                                                                                                                                   |
| wm.program.sleep               | スリープキーバインドを実行するコマンド                                                                                                                                                         |
| wm.program.settings            | コントロールパネルキーバインドを実行するコマンド                                                                                                                                                 |
| wm.program.display             | ディスプレイ設定キーバインドを実行するコマンド                                                                                                                                              |
| wm.program.wifi                | WiFi設定キーバインドを実行するコマンド                                                                                                                                                 |
| wm.program.bluetooth           | Bluetooth設定キーバインドを実行するコマンド                                                                                                                                            |
| wm.program.files               | ファイルブラウザキーバインドを実行するコマンド                                                                                                                                                  |
| wm.program.notification_ui     | 通知UIキーバインドを実行するコマンド                                                                                                                                               |
| wm.program.1...wm.program.3 | i3と開始するオプションのプログラム                                                                                                                                                                  |
| i3xrocks.media-player.player-args | `media-player`i3xrocksブロックにおけるオプションの引数。これらの引数は`playerctl`に直接与えられます。例えば、`-p spotify`という値を設定すると、バーがSpotifyを表示・操作するだけになります。 |
| i3xrocks.action.mic.left | マイクステータスインジケーターで左クリック時に実行するオプションのコマンド |
| i3xrocks.action.mic.middle | マイクステータスインジケーターで中クリック時に実行するオプションのコマンド |
| i3xrocks.action.mic.right | マイクステータスインジケーターで右クリック時に実行するオーバーライドコマンド |
| i3xrocks.action.mic.scrollup | マイクステータスインジケーターで上方向にスクロール時に実行するオーバーライドコマンド |
| i3xrocks.action.mic.scrolldn | マイクステータスインジケーターで下方向にスクロール時に実行するオーバーライドコマンド |
| i3xrocks.action.volume.left | 音量ステータスインジケーターで左クリック時に実行するオプションのコマンド |
| i3xrocks.action.volume.middle | 音量ステータスインジケーターで中クリック時に実行するオプションのコマンド |
| i3xrocks.action.volume.right | 音量ステータスインジケーターで右クリック時に実行するオーバーライドコマンド |
| i3xrocks.action.volume.scrollup | 音量ステータスインジケーターで下方向にスクロール時に実行するオーバーライドコマンド |
| i3xrocks.action.volume.scrolldn | 音量ステータスインジケーターで上方向にスクロール時に実行するオーバーライドコマンド |
| i3xrocks.action.disk-capacity | ディスク容量ステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.cpu-usage | CPU使用率ステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.memory | メモリステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.wifi | WiFiステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.weather | 天気ステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.media-player.left | メディアプレーヤーインジケーターを左クリック時に実行するオプションのコマンド |
| i3xrocks.action.media-player.middle | メディアプレーヤーインジケーターを中クリック時に実行するオプションのコマンド |
| i3xrocks.action.media-player.right | メディアプレーヤーインジケーターを右クリック時に実行するオーバーライドコマンド |
| i3xrocks.action.time | 日時ステータスインジケーターで実行するオーバーライドコマンド |
| i3xrocks.action.battery | バッテリーステータスインジケーターで実行するオーバーライドコマンド |
