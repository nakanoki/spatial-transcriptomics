# spatial-transcriptomics

10x Genomics Visium（ヒト乳がん・新鮮凍結）の空間トランスクリプトミクス解析。

## 依存関係

uv を使う。`uv run <cmd>` / `uv add <pkg>`。システム Python に install しない。

`pyyaml` / `leidenalg` / `squidpy` は `pyproject.toml` の `dependencies` に宣言済み。

空間データの読み込み・可視化は **squidpy を使う**（`sq.read.visium` / `sq.pl.spatial_scatter`）。
scanpy の `sc.read_visium` / `sc.pl.spatial` は squidpy へ移管され将来削除されるため、新規に使わない。

## データ

- `data/raw/` は読み取り専用。加工結果は `data/processed/` に出す。
- 合計 5.5GB。再取得コストが高いので、移動・削除・上書きをしない。
- **`data/` は .gitignore されており、リポジトリに含まれない。**
  GitHub Actions のランナー上には存在しない。

内訳（CI で使えるかの判断に使う）:

| ファイル | サイズ | 用途 |
| --- | ---: | --- |
| `Visium_Human_Breast_Cancer_image.tif` | 3.2G | H&E 画像。重ね合わせ表示用 |
| `..._cloupe.cloupe` | 1.8G | Loupe Browser 専用。Python 解析では使わない |
| `..._molecule_info.h5` | 200M | 用途次第 |
| `filtered_feature_bc_matrix/` ほか | 各 25〜60M | **解析の本体はここ** |

## コード構成

- `.ipynb` にロジックを書き溜めない。処理は `src/` に関数として置き、notebook からは import する。
- パラメータは `config/config.yaml` と `src/` の dataclass（`PreprocessParams` 等）で管理する。
  コードに直書きしない。
- notebook の出力セル（巨大な DataFrame、長いログ）をそのまま読み込まない。

## 完了の判断

テスト基盤が無い。動作確認をしていない変更を「完了」と報告しない。
最低限、変更したモジュールが import できること、対象スクリプトが実行できることを確認する。

CI 上ではデータが無いため実行確認ができない。その場合は完了扱いにせず、
「ローカルでの確認が必要」と明記する。
