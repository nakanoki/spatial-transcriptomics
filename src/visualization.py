"""
可視化ヘルパー（UMAP / QC / Spatial）
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import scanpy as sc


def plot_qc_violin(
    adata: sc.AnnData,
    keys: Iterable[str] = ("total_counts", "n_genes_by_counts", "pct_counts_mt"),
    groupby: str | None = None,
    save: str | Path | None = None,
) -> None:
    """
    QC 指標の violin plot。
    """
    sc.pl.violin(adata, keys=list(keys), groupby=groupby, multi_panel=True, save=None)
    if save is not None:
        # scanpy の save は figures 設定に依存するため、明示保存は呼び出し側で対応する想定
        raise NotImplementedError("明示的な保存は呼び出し側で実装してください（scanpy設定に依存）")


def plot_umap(
    adata: sc.AnnData,
    color: str | list[str] = "leiden",
    save: str | Path | None = None,
) -> None:
    """
    UMAP プロット（事前に sc.tl.umap が必要）。
    """
    sc.pl.umap(adata, color=color, save=None)
    if save is not None:
        raise NotImplementedError("明示的な保存は呼び出し側で実装してください（scanpy設定に依存）")


def plot_spatial(
    adata: sc.AnnData,
    color: str | list[str] = "leiden",
    library_id: str | None = None,
    img_key: str | None = None,
    alpha_img: float = 0.8,
    size: float = 1.2,
    save: str | Path | None = None,
) -> None:
    """
    Visium 形式の空間プロット（`sc.read_visium` を想定）。
    """
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] がありません（Visium 形式の AnnData を想定）")

    sc.pl.spatial(
        adata,
        color=color,
        library_id=library_id,
        img_key=img_key,
        alpha_img=alpha_img,
        size=size,
        save=None,
    )
    if save is not None:
        raise NotImplementedError("明示的な保存は呼び出し側で実装してください（scanpy設定に依存）")

