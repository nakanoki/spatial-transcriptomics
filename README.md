# spatial-transcriptomics

10x Genomics Visium（ヒト乳がん組織・Fresh Frozen）の空間トランスクリプトミクス解析プロジェクトです。

- **何を**: 乳がん組織切片の Visium 空間トランスクリプトームデータに対して、QC・正規化・次元削減・クラスタリング・マーカー遺伝子（DEG）探索・空間的な細胞間相互作用解析までを一貫して行うパイプラインを実装しています。
- **どのデータで**: 10x Genomics が公開している Human Breast Cancer（Visium, Fresh Frozen, Whole Transcriptome）データセットを使用します（詳細は [データソース](#データソース) を参照）。
- **何を明らかにするか**: 組織内の遺伝子発現に基づくクラスター（腫瘍領域・間質領域など）を同定し、それらの空間的な分布パターンやクラスター間の近接・相互作用を可視化・定量することを目指します。
- 実装は `scanpy` / `squidpy` をベースに、`config/config.yaml` による設定管理と `src/` 以下の再利用可能なモジュールで構成しています（ロジックを `.ipynb` に直書きしない方針、詳細は [CLAUDE.md](./CLAUDE.md) を参照）。

## ディレクトリ構成

```
config/           # 設定ファイル（パスやパラメータなど）
data/             # データ関連フォルダ（.gitignore 対象。下記「データソース」参照）
├─ processed/     # 前処理後のデータ
└─ raw/           # 生データ
    └─ breast_cancer/
        ├─ Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5   # Space Ranger 出力 H5 ファイル
        └─ spatial/                                                   # 空間座標・画像情報など
notebooks/        # 解析用のJupyterノートブック（ロジックは src/ から import して使う）
results/          # 解析結果（グラフ、表、レポートなど）
scripts/          # 実行用スクリプト（バッチ処理や再現用スクリプト）
src/              # 解析用のPythonコード（再利用可能な関数・dataclass パラメータ）
```

## データソース

- **名称**: Human Breast Cancer（Visium, Fresh Frozen, Whole Transcriptome）
- **URL**: https://www.10xgenomics.com/jp/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard
- **取得日**: `data/raw/` にダウンロードした日付をここに記入してください（このリポジトリでは実データを扱っていないため未記入）
- **ライセンス**: 10x Genomics が公開する Datasets の利用規約に準拠します。再配布・商用利用の可否を含む正確な条件は、上記 URL のページに記載された Terms of Use を参照してください。本リポジトリはデータそのものを再配布しません。
- **保管先**: S3（予定）。生データ・中間生成物・結果などの大きいファイルは S3 に置き、リポジトリには含めません。
- **備考**: Space Ranger 出力一式のうち、Python 解析で使うのは主に `filtered_feature_bc_matrix/`（または `.h5`）と `spatial/` です。H&E 画像（`.tif`）は重ね合わせ表示用、`.cloupe` は Loupe Browser 専用で Python 解析では使いません。詳細は [CLAUDE.md](./CLAUDE.md) のファイル一覧を参照してください。

### `data/` の取得手順

`data/` は `.gitignore` 対象のため、このリポジトリには含まれません（GitHub Actions のランナー上にも存在しません）。解析を再現するには、以下の手順でローカルにデータを配置してください。

1. 上記 URL から Space Ranger 出力一式（`filtered_feature_bc_matrix.h5` または `filtered_feature_bc_matrix/`、`spatial/` を含む）をダウンロードする。
2. `data/raw/<sample_id>/` （例: `data/raw/breast_cancer/`）以下に、ダウンロードしたファイル・フォルダをそのまま配置する。
3. `config/config.yaml` の `inputs.sample_id`（または `inputs.spaceranger_out_dir`）に配置先のディレクトリ名を設定する。

```yaml
inputs:
  sample_id: "breast_cancer"   # data/raw/breast_cancer を参照する場合
```

## セットアップ（uv）

依存管理には [uv](https://docs.astral.sh/uv/) を使用します。システム Python には install しません。

```bash
# 1. 依存関係のインストール
uv sync

# 2. データの配置（前節「data/ の取得手順」を参照）

# 3. パイプラインの実行
uv run python -c "
from pathlib import Path
from src.pipeline import run_pipeline
run_pipeline(Path('.'))
"
```

ノートブックで対話的に確認する場合:

```bash
uv run jupyter lab
# notebooks/visium_pipeline.ipynb を開く
```

パラメータ（QC 閾値・クラスタリング解像度など）は `config/config.yaml` と `src/` 各モジュールの dataclass（`PreprocessParams` / `DimReduceParams` / `ClusteringParams` 等）で管理しています。コードへの直書きはしません。

### 既知の制約

- `pyproject.toml` に宣言が漏れている依存関係があります（`pyyaml`, `leidenalg`）。`config.yaml` の `clustering.method: leiden` は現状 `leidenalg` が lock ファイルに含まれていないため失敗します。
- `squidpy`（空間近傍・cell-cell interaction 解析に使用）は任意導入で、未導入でもパイプライン本体は止まりません（`_require_squidpy()` で保護）。
- 詳細は [CLAUDE.md](./CLAUDE.md) を参照してください。

## `src/` モジュール一覧

| モジュール | 役割 |
| --- | --- |
| `src/preprocessing.py` | Visium データの読み込み（`sc.read_visium`）、QC メトリクス付与、セル/遺伝子フィルタ、正規化・log1p・高変動遺伝子（HVG）選択 |
| `src/clustering.py` | PCA → 近傍グラフ → UMAP による次元削減、Leiden/Louvain クラスタリング、Squidpy による空間近傍グラフ構築とモランI（空間自己相関）計算 |
| `src/deg.py` | `sc.tl.rank_genes_groups` をラップしたマーカー遺伝子・DEG（差次的発現）解析 |
| `src/visualization.py` | QC violin plot、UMAP、空間プロット（Visium 重ね合わせ表示）の可視化ヘルパー |
| `src/interaction.py` | Squidpy によるクラスター間の近傍出現頻度（neighborhood enrichment）・共起（co-occurrence）解析と可視化 |
| `src/pipeline.py` | `config/config.yaml` の設定に基づき、前処理〜次元削減〜クラスタリング〜（任意で）空間解析〜結果出力までを実行する統合パイプライン（`run_pipeline`） |
| `config/__init__.py` | `config/config.yaml` の読み込みとパス・パラメータ解決のユーティリティ（`load_config` / `get_paths` / `get_qc_config`） |

## Docker（R + Python + Jupyter）

このリポジトリには、R と Python の両方が使える Jupyter 環境用の `Dockerfile` が含まれます。`uv` によるローカル実行の代替手段です。

### ビルド

```bash
# プロジェクトルートで実行
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

## ライセンス

このリポジトリのコードは [MIT License](./LICENSE) のもとで公開しています。データセット自体のライセンスは [データソース](#データソース) 節を参照してください（コードのライセンスとは別です）。

## メモ

- `data/` と `results/` は `.dockerignore` でビルドコンテキストから除外しています（ビルド高速化のため）。
