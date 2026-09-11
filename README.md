# spatial-transcriptomics

[日本語](#日本語) | [English](#english)

---

## 日本語

10x Genomics Visium（ヒト乳がん組織・Fresh Frozen）の空間トランスクリプトミクス解析プロジェクトである。

- **何を**: 乳がん組織切片の Visium 空間トランスクリプトームデータに対して、QC・正規化・次元削減・クラスタリング・マーカー遺伝子（DEG）探索・空間的な細胞間相互作用解析までを一貫して行うパイプラインを実装している。
- **どのデータで**: 10x Genomics が公開している Human Breast Cancer（Visium, Fresh Frozen, Whole Transcriptome）データセットを使用する（詳細は [データソース](#データソース) を参照）。
- **何を明らかにするか**: 組織内の遺伝子発現に基づくクラスター（腫瘍領域・間質領域など）を同定し、それらの空間的な分布パターンやクラスター間の近接・相互作用を可視化・定量することを目指す。
- 実装は `scanpy` / `squidpy` をベースに、`config/config.yaml` による設定管理と `src/` 以下の再利用可能なモジュールで構成する（ロジックを `.ipynb` に直書きしない方針、詳細は [CLAUDE.md](./CLAUDE.md) を参照）。

### ディレクトリ構成

```
config/           # 設定ファイル（パスやパラメータなど）
data/             # データ関連フォルダ（.gitignore 対象。下記「データソース」参照）
├─ processed/     # 前処理後のデータ
└─ raw/           # 生データ
    └─ breast_cancer/
        ├─ Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5   # Space Ranger 出力 H5 ファイル
        └─ spatial/                                                   # 空間座標・画像情報など
notebooks/        # 解析用の Jupyter ノートブック（ロジックは src/ から import して使う）
results/          # 解析結果（グラフ、表、レポートなど）
scripts/          # 実行用スクリプト（バッチ処理や再現用スクリプト）
src/              # 解析用の Python コード（再利用可能な関数・dataclass パラメータ）
```

### データソース

- **名称**: Human Breast Cancer（Visium, Fresh Frozen, Whole Transcriptome）
- **URL**: https://www.10xgenomics.com/jp/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard
- **取得日**: `data/raw/` にダウンロードした日付をここに記入する（このリポジトリでは実データを扱っていないため未記入）
- **ライセンス**: 10x Genomics が公開する Datasets の利用規約に準拠する。再配布・商用利用の可否を含む正確な条件は、上記 URL のページに記載された Terms of Use を参照すること。本リポジトリはデータそのものを再配布しない。
- **保管先**: S3（予定）。生データ・中間生成物・結果などの大きいファイルは S3 に置き、リポジトリには含めない。
- **備考**: Space Ranger 出力一式のうち、Python 解析で使うのは主に `filtered_feature_bc_matrix/`（または `.h5`）と `spatial/` である。H&E 画像（`.tif`）は重ね合わせ表示用、`.cloupe` は Loupe Browser 専用で Python 解析では使わない。詳細は [CLAUDE.md](./CLAUDE.md) のファイル一覧を参照すること。

#### `data/` の取得手順

`data/` は `.gitignore` 対象のため、このリポジトリには含まれない（GitHub Actions のランナー上にも存在しない）。解析を再現するには、以下の手順でローカルにデータを配置する。

1. 上記 URL から Space Ranger 出力一式（`filtered_feature_bc_matrix.h5` または `filtered_feature_bc_matrix/`、`spatial/` を含む）をダウンロードする。
2. `data/raw/<sample_id>/` （例: `data/raw/breast_cancer/`）以下に、ダウンロードしたファイル・フォルダをそのまま配置する。
3. `config/config.yaml` の `inputs.sample_id`（または `inputs.spaceranger_out_dir`）に配置先のディレクトリ名を設定する。

```yaml
inputs:
  sample_id: "breast_cancer"   # data/raw/breast_cancer を参照する場合
```

### セットアップ（uv）

依存管理には [uv](https://docs.astral.sh/uv/) を使用する。システム Python には install しない。

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

パラメータ（QC 閾値・クラスタリング解像度など）は `config/config.yaml` と `src/` 各モジュールの dataclass（`PreprocessParams` / `DimReduceParams` / `ClusteringParams` 等）で管理する。コードへの直書きはしない。

#### 既知の制約

- `pyproject.toml` に宣言が漏れている依存関係がある（`pyyaml`, `leidenalg`）。`config.yaml` の `clustering.method: leiden` は現状 `leidenalg` が lock ファイルに含まれていないため失敗する。
- `squidpy`（空間近傍・cell-cell interaction 解析に使用）は任意導入で、未導入でもパイプライン本体は止まらない（`_require_squidpy()` で保護）。
- 詳細は [CLAUDE.md](./CLAUDE.md) を参照すること。

### `src/` モジュール一覧

| モジュール | 役割 |
| --- | --- |
| `src/preprocessing.py` | Visium データの読み込み（`sc.read_visium`）、QC メトリクス付与、セル/遺伝子フィルタ、正規化・log1p・高変動遺伝子（HVG）選択 |
| `src/clustering.py` | PCA → 近傍グラフ → UMAP による次元削減、Leiden/Louvain クラスタリング、Squidpy による空間近傍グラフ構築とモランI（空間自己相関）計算 |
| `src/deg.py` | `sc.tl.rank_genes_groups` をラップしたマーカー遺伝子・DEG（差次的発現）解析 |
| `src/visualization.py` | QC violin plot、UMAP、空間プロット（Visium 重ね合わせ表示）の可視化ヘルパー |
| `src/interaction.py` | Squidpy によるクラスター間の近傍出現頻度（neighborhood enrichment）・共起（co-occurrence）解析と可視化 |
| `src/pipeline.py` | `config/config.yaml` の設定に基づき、前処理〜次元削減〜クラスタリング〜（任意で）空間解析〜結果出力までを実行する統合パイプライン（`run_pipeline`） |
| `config/__init__.py` | `config/config.yaml` の読み込みとパス・パラメータ解決のユーティリティ（`load_config` / `get_paths` / `get_qc_config`） |

### Docker（R + Python + Jupyter）

このリポジトリには、R と Python の両方が使える Jupyter 環境用の `Dockerfile` が含まれる。`uv` によるローカル実行の代替手段である。

#### ビルド

```bash
# プロジェクトルートで実行
docker build -t spatial-rpy .
```

#### 起動（ノートブック実行）

```bash
docker run --rm -p 8888:8888 -v "$PWD":/work spatial-rpy
```

起動ログに表示される **token 付き URL** をブラウザで開く。

### 解析環境（EC2）

（未確定。決まり次第ここに追記する）

- **リージョン**:
- **AMI**:
- **インスタンスタイプ**:
- **ストレージ**:（例: EBS の容量/種類）
- **ネットワーク**:（例: VPC / セキュリティグループ）
- **データ配置**:（例: S3 バケット名・prefix、ローカルに同期するディレクトリ）
- **実行方法**:（例: Docker を使う / 直接環境構築する）
- **備考**:

### ライセンス

このリポジトリのコードは [MIT License](./LICENSE) のもとで公開している。データセット自体のライセンスは [データソース](#データソース) 節を参照すること（コードのライセンスとは別である）。

### メモ

- `data/` と `results/` は `.dockerignore` でビルドコンテキストから除外している（ビルド高速化のため）。

---

## English

A spatial transcriptomics analysis project using 10x Genomics Visium data (human breast cancer tissue, fresh frozen).

- **What**: Implements an end-to-end pipeline over Visium spatial transcriptome data from breast cancer tissue sections — QC, normalization, dimensionality reduction, clustering, marker gene (DEG) discovery, and spatial cell-cell interaction analysis.
- **What data**: Uses the publicly available Human Breast Cancer (Visium, Fresh Frozen, Whole Transcriptome) dataset from 10x Genomics (see [Data Source](#data-source) for details).
- **Goal**: Identify expression-based clusters within the tissue (e.g. tumor regions, stromal regions), and visualize/quantify their spatial distribution patterns and inter-cluster proximity/interactions.
- Built on `scanpy` / `squidpy`, with configuration managed via `config/config.yaml` and reusable modules under `src/` (logic is not written directly in `.ipynb` files — see [CLAUDE.md](./CLAUDE.md) for details).

### Directory Structure

```
config/           # Config files (paths, parameters, etc.)
data/             # Data folder (gitignored; see "Data Source" below)
├─ processed/     # Preprocessed data
└─ raw/           # Raw data
    └─ breast_cancer/
        ├─ Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5   # Space Ranger output H5 file
        └─ spatial/                                                   # Spatial coordinates / image info, etc.
notebooks/        # Jupyter notebooks for analysis (imports logic from src/)
results/          # Analysis outputs (plots, tables, reports, etc.)
scripts/          # Executable scripts (batch jobs, reproduction scripts)
src/              # Python code for analysis (reusable functions / dataclass parameters)
```

### Data Source

- **Name**: Human Breast Cancer (Visium, Fresh Frozen, Whole Transcriptome)
- **URL**: https://www.10xgenomics.com/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard
- **Retrieved on**: fill in the date the data was downloaded into `data/raw/` here (left blank in this repository since no real data is handled here)
- **License**: Governed by the Terms of Use published on 10x Genomics' Datasets page. Refer to the Terms of Use at the URL above for the exact conditions, including redistribution and commercial use. This repository does not redistribute the data itself.
- **Storage**: S3 (planned). Large files such as raw data, intermediate outputs, and results are stored on S3 and are not included in this repository.
- **Notes**: Of the full Space Ranger output, Python analysis mainly uses `filtered_feature_bc_matrix/` (or `.h5`) and `spatial/`. The H&E image (`.tif`) is used for overlay visualization, and `.cloupe` is for Loupe Browser only — it is not used in Python analysis. See the file list in [CLAUDE.md](./CLAUDE.md) for details.

#### Getting `data/`

`data/` is gitignored and is not included in this repository (nor does it exist on GitHub Actions runners). To reproduce the analysis, place the data locally following these steps:

1. Download the full Space Ranger output (`filtered_feature_bc_matrix.h5` or `filtered_feature_bc_matrix/`, including `spatial/`) from the URL above.
2. Place the downloaded files/folders as-is under `data/raw/<sample_id>/` (e.g. `data/raw/breast_cancer/`).
3. Set `inputs.sample_id` (or `inputs.spaceranger_out_dir`) in `config/config.yaml` to the destination directory name.

```yaml
inputs:
  sample_id: "breast_cancer"   # when referencing data/raw/breast_cancer
```

### Setup (uv)

Dependency management uses [uv](https://docs.astral.sh/uv/). Nothing is installed into the system Python.

```bash
# 1. Install dependencies
uv sync

# 2. Place the data (see "Getting data/" above)

# 3. Run the pipeline
uv run python -c "
from pathlib import Path
from src.pipeline import run_pipeline
run_pipeline(Path('.'))
"
```

To explore interactively via notebook:

```bash
uv run jupyter lab
# open notebooks/visium_pipeline.ipynb
```

Parameters (QC thresholds, clustering resolution, etc.) are managed via `config/config.yaml` and dataclasses in each `src/` module (`PreprocessParams` / `DimReduceParams` / `ClusteringParams`, etc.), not hardcoded.

#### Known Limitations

- `pyproject.toml` is missing some dependency declarations (`pyyaml`, `leidenalg`). Since `leidenalg` is absent from the lock file, `clustering.method: leiden` in `config.yaml` currently fails.
- `squidpy` (used for spatial neighborhood / cell-cell interaction analysis) is an optional dependency; the core pipeline still runs without it (guarded by `_require_squidpy()`).
- See [CLAUDE.md](./CLAUDE.md) for details.

### `src/` Module Overview

| Module | Role |
| --- | --- |
| `src/preprocessing.py` | Loads Visium data (`sc.read_visium`), attaches QC metrics, filters cells/genes, normalizes, log1p, and selects highly variable genes (HVGs) |
| `src/clustering.py` | Dimensionality reduction via PCA → neighbor graph → UMAP, Leiden/Louvain clustering, spatial neighbor graph construction via Squidpy, and Moran's I (spatial autocorrelation) computation |
| `src/deg.py` | Marker gene / DEG (differential expression) analysis wrapping `sc.tl.rank_genes_groups` |
| `src/visualization.py` | Visualization helpers for QC violin plots, UMAP, and spatial plots (Visium overlay) |
| `src/interaction.py` | Neighborhood enrichment and co-occurrence analysis/visualization between clusters via Squidpy |
| `src/pipeline.py` | Integrated pipeline (`run_pipeline`) that runs preprocessing → dimensionality reduction → clustering → (optionally) spatial analysis → output, based on `config/config.yaml` |
| `config/__init__.py` | Utilities for loading `config/config.yaml` and resolving paths/parameters (`load_config` / `get_paths` / `get_qc_config`) |

### Docker (R + Python + Jupyter)

This repository includes a `Dockerfile` for a Jupyter environment with both R and Python available, as an alternative to running locally with `uv`.

#### Build

```bash
# run from the project root
docker build -t spatial-rpy .
```

#### Run (launch notebook)

```bash
docker run --rm -p 8888:8888 -v "$PWD":/work spatial-rpy
```

Open the **URL with the token** shown in the startup log in a browser.

### Analysis Environment (EC2)

(Not yet decided; will be filled in once determined)

- **Region**:
- **AMI**:
- **Instance type**:
- **Storage**: (e.g. EBS size/type)
- **Network**: (e.g. VPC / security group)
- **Data placement**: (e.g. S3 bucket/prefix, local sync directory)
- **How to run**: (e.g. Docker vs. direct environment setup)
- **Notes**:

### License

The code in this repository is released under the [MIT License](./LICENSE). See the [Data Source](#data-source) section for the dataset's own license (separate from the code license).

### Notes

- `data/` and `results/` are excluded from the Docker build context via `.dockerignore` (to speed up builds).
