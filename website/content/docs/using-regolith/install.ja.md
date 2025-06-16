---
title: "インストール"
linkTitle: "インストール"
weight: 1
description: >
  Regolithのインストール
prev: /docs/using-regolith
---

## Upgrading to 3.x

以前のバージョンのRegolithからアップグレードする場合は[リリースノート](/docs/reference/Releases/regolith-3.2-release-notes) を、特に2.xからの場合は[移行ガイド](/docs/reference/Releases/regolith-3.0-release-notes/#migration-guide)を参照してください。

## サポートするOS

以下は、Regolithデスクトップのインストールを取り扱っている特定のLinuxベースのOSごとに分割されています。

### Ubuntu

{{< tabs items="Ubuntu 24.04,Ubuntu 22.04" >}}
{{< tab >}}
Regolithはシステムパッケージとしてインストールすることができます。これは簡単に更新や削除ができ、一貫することができます。
下記のステップでは、Regolithのパッケージリポジトリからパッケージを読み込み可能にして、デスクトップパッケージをインストールできるようにシステムを設定する方法について記しています。

1. ローカルの`apt`にRegolithの公開鍵を登録する。

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. ローカルの`apt`にリポジトリのURLを追加する。

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/ubuntu/stable noble v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. `apt`を更新してRegolithをインストールする。

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. システムを再起動する。

ログインマネージャーは新しいデスクトップセッションを識別するため、再起動が必要になります。システムの再起動はログインマネージャーを再起動するための簡単な方法です。

{{< callout type="info" >}}
ARMベースのシステムにインストールするときは、上記の内容から2箇所の`amd64`を`arm64`に置き換えてください。
{{< /callout >}}
{{< /tab >}}

{{< tab >}}
Regolithはシステムパッケージとしてインストールすることができます。これは簡単に更新や削除ができ、一貫することができます。
下記のステップでは、Regolithのパッケージリポジトリからパッケージを読み込み可能にして、デスクトップパッケージをインストールできるようにシステムを設定する方法について記しています([こちらでダウンロード可能なスクリプト](/install-release-ubuntu-22.04-amd64.txt)が利用できます)。

1. ローカルの`apt`にRegolithの公開鍵を登録する。

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. ローカルの`apt`にリポジトリのURLを追加する。

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/ubuntu/stable jammy v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. `apt`を更新してRegolithをインストールする。

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. システムを再起動する。

ログインマネージャーは新しいデスクトップセッションを識別するため、再起動が必要になります。システムの再起動はログインマネージャーを再起動するための簡単な方法です。

{{< callout type="info" >}}
ARMベースのシステムにインストールするときは、上記の内容から2箇所の`amd64`を`arm64`に置き換えてください。
{{< /callout >}}
{{< /tab >}}
{{< /tabs >}}

### Debian

{{< tabs items="Debian Testing,Debian Bookworm" >}}
{{< tab >}}
"Testing"リリースのDebianユーザーであれば、Regolithをインストールすることができますが、Debian Testingの性質上"release"バージョンは存在しません。Debian TestingインスタンスにRegolithをインストールするには、`apt`コマンドのワンライナーの`release`の代わりに`testing`を含める必要があります。例：

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/debian/testing testing main" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

それ以外のインストールの手順はDebian公式リリースと同じです。壊れた依存関係やコンポーネントと時々遭遇することに注意してください。安定性が必要である場合は、DebianやUbuntuの公式リリースにインストールすることを検討してください。
{{< /tab >}}

{{< tab >}}
Regolithはシステムパッケージとしてインストールすることができます。これは簡単に更新や削除ができ、一貫することができます。
下記のステップでは、Regolithのパッケージリポジトリからパッケージを読み込み可能にして、デスクトップパッケージをインストールできるようにシステムを設定する方法について記しています([こちらでダウンロード可能なスクリプト](/install-release-debian-12-amd64.txt)が利用できます)。

1. ローカルの`apt`にRegolithの公開鍵を登録する。

   ```bash
   wget -qO - https://archive.regolith-desktop.com/regolith.key | \
   gpg --dearmor | sudo tee /usr/share/keyrings/regolith-archive-keyring.gpg > /dev/null
   ```

1. ローカルの`apt`にリポジトリのURLを追加する。

   ```bash
   echo deb "[arch=amd64 signed-by=/usr/share/keyrings/regolith-archive-keyring.gpg] \
   https://archive.regolith-desktop.com/debian/stable bookworm v3.2" | \
   sudo tee /etc/apt/sources.list.d/regolith.list
   ```

1. `apt`を更新してRegolithをインストールする。

   ```bash
   sudo apt update
   sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
   ```
1. システムを再起動する。

ログインマネージャーは新しいデスクトップセッションを識別するため、再起動が必要になります。システムの再起動はログインマネージャーを再起動するための簡単な方法です。

{{< callout type="info" >}}
ARMベースのシステムにインストールするときは、上記の内容から2箇所の`amd64`を`arm64`に置き換えてください。
{{< /callout >}}

{{< callout type="info" >}}
`regolith-desktop`パッケージは、多くのハードウェアで動作する必要最小限の設定でインストールします。追加で高級なパッケージやもっと便利な環境を追加するには、[新規ユーザーに推奨されるパッケージ]({{< ref "docs/using-regolith/configuration.md#recommended-packages-for-new-users" >}})を参照してください。
{{< /callout >}}

{{< callout type="info" >}}
`regolith-compositor-picom-glx`コンポジターは多くのコンピューターで動作します。もしドライバや視覚的な問題に遭遇した場合は、[他のコンポジターをお試しください]({{< ref "docs/howtos/customize-compositor" >}})。
{{< /callout >}}
{{< /tab >}}
{{< /tabs >}}

## インストール時に選択するパッケージ

Regolith 3.0では、複数のセッション(X11: `regolith-session-flashback`, Wayland: `regolith-session-sway`)から選択することができるようになりました。このため、Regolith 3.0をインストールするときは、使用するセッションを一方か両方を指定します。さらに、Regolith 3.0からはインストール時にRegolithの外観を指定できるようになったことで、後で設定する必要があるセットアップが節約されます。DebianベースのシステムのX11で推奨される基本的なインストールは以下の通りです。

```bash
sudo apt install regolith-desktop regolith-session-flashback regolith-look-lascaille
#                ^-- ﾍﾞｰｽ ﾊﾟｯｹｰｼﾞ   ^-- ｾｯｼｮﾝ                   ^ -- 外観 
```

* セッションについて、次の内の一方か両方を選ぶことができます：`regolith-session-flashback` (X11), `regolith-session-sway` (Wayland)
* 外観について、次の外観パッケージから*1つ*を選ぶことができます。
```text
regolith-look-ayu-dark
regolith-look-ayu-mirage
regolith-look-ayu
regolith-look-blackhole
regolith-look-default-loader
regolith-look-default
regolith-look-dracula
regolith-look-gruvbox
regolith-look-i3-default
regolith-look-lascaille
regolith-look-nevil
regolith-look-nord
regolith-look-solarized-dark
```

上記の推奨事項の例と代替案として、ここでは、Swayセッションに外観"Nord"を適用してインストールする`apt`コマンドを紹介します。

```bash
sudo apt install regolith-desktop regolith-session-sway regolith-look-nord
```

## Anatomy of the Regolith `apt` Config Line

```text
deb [arch=amd64] https://archive.regolith-desktop.com/ubuntu/unstable noble main
          |                                           |      |        |     |
          |                                           |      |        |     * Repo Component
          |                                           |      |        * ディストロバージョン/コードネーム
          |                                           |      * Repo Suite
          |                                           * ディストロ名
          * アーキテクチャ
```

### アーキテクチャ

| 名前 | 概要 |
|-----|:-----|
| `amd64` | 64-bit x86アーキテクチャ |
| `arm64` | 64-bit Armアーキテクチャ |

### ディストロ名

| 名前 | 概要 |
|-----|:-----|
| `ubuntu` | Ubuntu Linuxディストリビューション |
| `debian` | Debian Linuxディストリビューション |

### Repo Suites and Components

| Regolith 段階 | Repo Suite | Repo Component | 概要 |
|:-------------|------------|:--------------:|:-----|
| 実験的 | `experimental` | `main` | 実験的内容をテストする段階 |
| 不安定 | `unstable` | `main` | 開発内容の初期テスト段階 |
| テスト | `testing` | `main` | 開発内容の後期テスト段階 |
| リリース X.Y | `stable` | `vX.Y` | 特定の公式リリース段階 |
| 現在のリリース | `stable` | `main` | 最新のリリース段階（リリース間で変動します） |

### ディストロバージョン/コードネーム

これらのラベルは各上流コミュニティによって決定されます。例：`noble`, `bookworm`, `focal`


### 例

|概要 | APT Line |
|:---|:---------|
| `amd64`のUbuntu 22.04で3.2リリースのRegolithを使用する | `deb [arch=amd64] https://archive.regolith-desktop.com/ubuntu/stable jammy v3.2` |
| `amd64`のDebian 12で最新リリースのRegolithを使用する | `deb [arch=amd64] https://archive.regolith-desktop.com/debian/stable bookworm main` |
| `arm64`のDebian TestingでRegolithのテストリポジトリを使用する | `deb [arch=arm64] https://archive.regolith-desktop.com/debian/testing testing main` |


## リリースの変更ポリシー

 From Regolith 3.0 to 3.2 (inclusive) Regolith 3.2以降では、すべてのリリースにおいて、パッケージリポジトリのURLに一意の名前が使用されます。これは、ユーザーが、新しいリリースにアップグレードするタイミングを完全に制御できることを意味します。パッケージマネージャーに常に最新版をインストールしてほしいユーザーには、`release-current` という特別な段階が提供されています。

From Regolith 3.3 onward, all releases will be published in the common archive repository separated by a high level distro (e.g. `debian`, `ubuntu`) folder and corresponding codename subfolders. The structure is split furthermore into different suites and components. As such, users wishing to have their package manager always install the latest version, they should use `main` component of `stable` suite. Otherwise a fixed version component (for example `v3.3`, `v3.2`, etc) can be used alongside `stable` component.

Note that `experimental`, `unstable`, and `testing` suites only have `main` component.
