# spatial-transcriptomics

Spatial transcriptomics 解析用プロジェクト。

## ディレクトリ構成

- `notebooks/`: 解析ノートブック
- `src/`: 解析用コード
- `data/raw/`: 生データ（大きいので Git 管理しない想定）
- `data/processed/`: 前処理後データ
- `results/`: 解析結果
- `scripts/`: 実行スクリプト
- `config/`: 設定ファイル

## データソース

https://www.10xgenomics.com/jp/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard


- **保管先**: S3（予定）
  - **方針**: 生データ・中間生成物・結果などの大きいファイルは S3 に置き、リポジトリには含めない
  - **ローカル配置**: 必要に応じて `data/raw/` / `data/processed/` にダウンロードして利用

- **名称**:
  - **URL**:
  - **取得日**:
  - **ライセンス**:
  - **備考**:

## Docker（R + Python + Jupyter）

このリポジトリには、R と Python の両方が使える Jupyter 環境用の `Dockerfile` が含まれます。

### ビルド

```bash
cd /Users/kie/projects/spatial-transcriptomics
docker build -t spatial-rpy .
```

### 起動（ノートブック実行）

```bash
docker run --rm -p 8888:8888 -v "$PWD":/work spatial-rpy
```

起動ログに表示される **token 付きURL** をブラウザで開いてください。

## 解析環境（EC2）

（未確定。決まり次第ここに追記）

- **リージョン**:
- **AMI**:
- **インスタンスタイプ**:
- **ストレージ**:（例: EBS の容量/種類）
- **ネットワーク**:（例: VPC / セキュリティグループ）
- **データ配置**:（例: S3 バケット名・prefix、ローカルに同期するディレクトリ）
- **実行方法**:（例: Docker を使う / 直接環境構築する）
- **備考**:

## メモ

- `data/` と `results/` は `.dockerignore` でビルドコンテキストから除外しています（ビルド高速化のため）。

