"""
可視化ヘルパー（UMAP / QC / Spatial）
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import scanpy as sc
import squidpy as sq


def plot_qc_violin(
    adata: sc.AnnData,
    keys: Iterable[str] = ("total_counts", "n_genes_by_counts", "pct_counts_mt"),
    groupby: str | None = None,
    save: str | Path | None = None,
    show: bool = False,
    dpi: int = 150,
):
    """
    QC 指標の violin plot。

    `save` にパスを渡すと、その場所に画像として保存する。
    notebook から対話的に見たい場合は `show=True`。
    """
    ret = sc.pl.violin(
        adata, keys=list(keys), groupby=groupby, multi_panel=True, show=show
    )
    if save is not None:
        import matplotlib.pyplot as plt

        save = Path(save)
        save.parent.mkdir(parents=True, exist_ok=True)
        # multi_panel=True は seaborn の FacetGrid(.fig)、それ以外は Axes(.figure) を返す
        fig = getattr(ret, "fig", None) or getattr(ret, "figure", None) or plt.gcf()
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return ret


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
    dpi: int = 150,
) -> None:
    """
    Visium 形式の空間プロット（`sq.read.visium` で読んだ AnnData を想定）。

    `save` にパスを渡すと、その場所に画像として保存する。
    """
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] がありません（Visium 形式の AnnData を想定）")

    sq.pl.spatial_scatter(
        adata,
        color=color,
        library_id=library_id,
        img_res_key=img_key,
        alpha_img=alpha_img,
        size=size,
    )
    if save is not None:
        import matplotlib.pyplot as plt

        save = Path(save)
        save.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        plt.close("all")

