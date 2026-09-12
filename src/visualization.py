"""
可視化ヘルパー（UMAP / QC / Spatial）
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import scanpy as sc
import squidpy as sq


def _set_title(title: str | None) -> None:
    """
    図全体のタイトルを設定する。`bbox_inches="tight"` での保存を前提に余白を取る。

    タイトルは英語で書くこと。matplotlib の既定フォントは日本語を持たず、
    豆腐（□）になる。日本語を出すには環境ごとにフォント設定が必要で、
    CI や他マシンで同じ図が再現できなくなる。
    サブプロット名と重なるため1行に収めること。
    """
    if not title:
        return
    import matplotlib.pyplot as plt

    plt.gcf().suptitle(title, fontsize=13, y=1.04)


def plot_qc_violin(
    adata: sc.AnnData,
    keys: Iterable[str] = ("total_counts", "n_genes_by_counts", "pct_counts_mt"),
    groupby: str | None = None,
    save: str | Path | None = None,
    show: bool = False,
    dpi: int = 150,
    title: str | None = None,
    ylims: dict[str, tuple[float, float]] | None = None,
):
    """
    QC 指標の violin plot。

    `save` にパスを渡すと、その場所に画像として保存する。
    notebook から対話的に見たい場合は `show=True`。

    `ylims` に `{指標名: (下限, 上限)}` を渡すと、その範囲に軸を固定する。
    フィルタ前後の図を並べて比較するときは、前者の範囲を両方に与えること。
    自動スケールのままでは軸が変わり、分布の変化を見誤る。
    """
    ret = sc.pl.violin(
        adata, keys=list(keys), groupby=groupby, multi_panel=True, show=show
    )
    if ylims:
        import matplotlib.pyplot as plt

        # multi_panel=True の軸は keys の順に並ぶ
        for ax, key in zip(plt.gcf().axes, list(keys)):
            if key in ylims:
                ax.set_ylim(*ylims[key])

    _set_title(title)
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
    img_alpha: float | None = None,
    size: float | None = None,
    save: str | Path | None = None,
    dpi: int = 150,
    crop: bool = True,
    crop_margin: float = 0.02,
    title: str | None = None,
) -> None:
    """
    Visium 形式の空間プロット（`sq.read.visium` で読んだ AnnData を想定）。

    `save` にパスを渡すと、その場所に画像として保存する。

    squidpy は既定でスライド画像の全体を描画するため、`crop=True` のときは
    スポットの存在範囲（+ `crop_margin` の余白）に切り詰める。
    組織以外の余白や、同一スライド上の別の組織片が写り込むのを防ぐ。
    """
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] がありません（Visium 形式の AnnData を想定）")

    # None の引数は渡さず squidpy の既定に委ねる。
    # 特に img_res_key は既定が "hires" で、None を明示すると画像とスポットの
    # 座標系がずれてスポットが描画されなくなる。
    kwargs: dict = {"color": color, "library_id": library_id}
    if img_key is not None:
        kwargs["img_res_key"] = img_key
    if img_alpha is not None:
        kwargs["img_alpha"] = img_alpha
    if size is not None:
        kwargs["size"] = size

    sq.pl.spatial_scatter(adata, **kwargs)
    if crop:
        _crop_to_spots(adata, library_id=library_id, img_key=img_key, margin=crop_margin)

    _set_title(title)

    if save is not None:
        import matplotlib.pyplot as plt

        save = Path(save)
        save.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save, dpi=dpi, bbox_inches="tight")
        plt.close("all")


def _crop_to_spots(
    adata: sc.AnnData,
    library_id: str | None = None,
    img_key: str | None = None,
    margin: float = 0.02,
) -> None:
    """
    描画済みの図の表示範囲を、スポットの存在範囲（+ margin）に詰める。

    squidpy は既定でスライド画像の全体を描画するため、組織以外の余白や
    同一スライド上の別の組織片まで写り込む。

    軸の座標系は画像のピクセル座標なので、`adata.obsm["spatial"]`
    （フル解像度座標）に scalefactor を掛けた値が、そのまま表示範囲になる。
    """
    if "spatial" not in adata.uns:
        return

    import matplotlib.pyplot as plt

    lib = library_id or next(iter(adata.uns["spatial"]))
    scalefactors = adata.uns["spatial"][lib].get("scalefactors", {})
    res = img_key or "hires"
    scale = scalefactors.get(f"tissue_{res}_scalef", 1.0)

    xy = adata.obsm["spatial"] * scale
    x0, y0 = xy.min(axis=0)
    x1, y1 = xy.max(axis=0)
    mx = (x1 - x0) * margin
    my = (y1 - y0) * margin

    for ax in plt.gcf().axes:
        if not ax.images:          # 画像を持つ軸だけが対象（カラーバーは除く）
            continue
        y_lo, y_hi = ax.get_ylim()
        ax.set_xlim(x0 - mx, x1 + mx)
        # 画像の軸は y が反転しているため、現在の向きを保って設定する
        if y_lo > y_hi:
            ax.set_ylim(y1 + my, y0 - my)
        else:
            ax.set_ylim(y0 - my, y1 + my)
