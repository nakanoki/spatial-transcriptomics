"""
Issue #2: 品質管理（QC）の実行スクリプト。

config/config.yaml の閾値で QC を行い、以下を出力する。

    results/qc/qc_violin_before.png   フィルタ前の QC 分布
    results/qc/qc_violin_after.png    フィルタ後の QC 分布
    results/qc/qc_spatial.png         QC 指標の空間分布
    data/processed/<sample_id>_qc.h5ad  フィルタ済みデータ

使い方:  uv run python scripts/run_qc.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 画面を持たない環境で実行するため
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import get_paths, load_config  # noqa: E402
from src.preprocessing import (  # noqa: E402
    add_qc_metrics,
    filter_cells_and_genes,
    read_visium_sample,
)
from src.visualization import plot_qc_violin, plot_spatial  # noqa: E402

QC_KEYS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]


def describe(adata: sc.AnnData, label: str) -> None:
    """QC 指標の分布を数値で出す（閾値を決める材料）。"""
    print(f"\n--- {label}: {adata.n_obs} spots x {adata.n_vars} genes ---")
    print(adata.obs[QC_KEYS].describe().round(2).to_string())


def main() -> None:
    cfg = load_config(PROJECT_ROOT)
    paths = get_paths(PROJECT_ROOT, cfg=cfg)
    qc = cfg["qc"]

    sample_id = cfg["inputs"]["sample_id"]
    sample_path = paths["data_raw"] / sample_id
    out_dir = paths["results"] / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = read_visium_sample(sample_path)
    add_qc_metrics(adata)
    describe(adata, "フィルタ前")
    plot_qc_violin(adata, keys=QC_KEYS, save=out_dir / "qc_violin_before.png")

    # QC 指標の空間分布（組織の端や剥離の影響を確認する）
    plot_spatial(
        adata,
        color=["total_counts", "n_genes_by_counts"],
        save=out_dir / "qc_spatial.png",
    )

    print(
        f"\n閾値: min_genes(per spot)={qc['min_genes_per_spot']}, "
        f"min_cells(per gene)={qc['min_cells']}, "
        f"pct_counts_mt_max={qc['pct_counts_mt_max']}"
    )
    adata = filter_cells_and_genes(
        adata,
        min_genes=qc["min_genes_per_spot"],
        min_cells=qc["min_cells"],
        pct_counts_mt_max=qc["pct_counts_mt_max"],
    )
    describe(adata, "フィルタ後")
    plot_qc_violin(adata, keys=QC_KEYS, save=out_dir / "qc_violin_after.png")

    processed = paths["data_processed"]
    processed.mkdir(parents=True, exist_ok=True)
    out_h5ad = processed / f"{sample_id}_qc.h5ad"
    adata.write_h5ad(out_h5ad)

    print(f"\n保存:\n  {out_h5ad}\n  {out_dir}/qc_violin_before.png"
          f"\n  {out_dir}/qc_violin_after.png\n  {out_dir}/qc_spatial.png")


if __name__ == "__main__":
    main()
