"""
Issue #3: 正規化・対数変換・高変動遺伝子（HVG）選択の実行スクリプト。

QC（scripts/run_qc.py）の出力を入力とし、以下を出力する。

    results/preprocess/normalization_effect.png   正規化前後のスポット深度
    results/preprocess/hvg_selection.png          HVG の選択結果
    data/processed/<sample_id>_preprocessed.h5ad   前処理済みデータ

使い方:  uv run python scripts/run_preprocess.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 画面を持たない環境で実行するため
import scanpy as sc  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import get_paths, load_config  # noqa: E402
from src.preprocessing import load_processed_h5ad, normalize_and_hvg  # noqa: E402
from src.visualization import plot_hvg, plot_normalization_effect  # noqa: E402


def main() -> None:
    cfg = load_config(PROJECT_ROOT)
    paths = get_paths(PROJECT_ROOT, cfg=cfg)

    sample_id = cfg["inputs"]["sample_id"]
    target_sum = float(cfg["normalization"]["target_sum"])
    n_top_genes = int(cfg["feature_selection"]["n_top_genes"])

    out_dir = paths["results"] / "preprocess"
    out_dir.mkdir(parents=True, exist_ok=True)

    qc_h5ad = paths["data_processed"] / f"{sample_id}_qc.h5ad"
    if not qc_h5ad.exists():
        raise FileNotFoundError(
            f"QC の出力が見つかりません: {qc_h5ad}\n"
            "先に `uv run python scripts/run_qc.py` を実行してください。"
        )
    adata = load_processed_h5ad(qc_h5ad)
    print(f"入力 (QC 済み): {adata.n_obs:,} spots x {adata.n_vars:,} genes")

    # 正規化の効果は、正規化前のカウントを持つ状態でのみ示せる
    plot_normalization_effect(
        adata,
        target_sum=target_sum,
        save=out_dir / "normalization_effect.png",
        title=(
            "Normalization equalizes per-spot sequencing depth  |  "
            f"target_sum={target_sum:g}"
        ),
    )

    adata = normalize_and_hvg(adata, target_sum=target_sum, n_top_genes=n_top_genes)
    print(f"正規化・HVG 選択後: {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    print(f"adata.raw（HVG subset 前）: {adata.raw.n_vars:,} genes")

    # HVG のプロットには subset 前の全遺伝子が必要
    plot_hvg(
        adata.raw.to_adata(),
        save=out_dir / "hvg_selection.png",
        title=f"Highly variable gene selection  |  top {n_top_genes:,} of {adata.raw.n_vars:,}",
    )

    out_h5ad = paths["data_processed"] / f"{sample_id}_preprocessed.h5ad"
    adata.write_h5ad(out_h5ad)

    print(
        f"\n保存:\n  {out_h5ad}"
        f"\n  {out_dir}/normalization_effect.png"
        f"\n  {out_dir}/hvg_selection.png"
    )


if __name__ == "__main__":
    main()
