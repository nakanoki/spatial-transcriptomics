"""
cell-cell interaction / 近接関係解析のための Squidpy ラッパー。
"""

from __future__ import annotations

from dataclasses import dataclass

import scanpy as sc
import squidpy as sq


@dataclass(frozen=True)
class InteractionParams:
    cluster_key: str = "clusters"
    copy: bool = False


def neighborhood_enrichment(
    adata: sc.AnnData,
    params: InteractionParams = InteractionParams(),
) -> sc.AnnData | None:
    """
    クラスター間の近傍出現頻度（Neighborhood enrichment）を計算する。

    結果は `adata.uns['nhood_enrichment']` に格納される。
    """
    return sq.gr.nhood_enrichment(
        adata,
        cluster_key=params.cluster_key,
        copy=params.copy,
    )


def co_occurrence(
    adata: sc.AnnData,
    params: InteractionParams = InteractionParams(),
) -> sc.AnnData | None:
    """
    距離スケールごとのクラスター共起（Co-occurrence）を計算する。

    結果は `adata.uns['co_occurrence']` に格納される。
    """
    return sq.gr.co_occurrence(
        adata,
        cluster_key=params.cluster_key,
        copy=params.copy,
    )


def plot_neighborhood_enrichment(
    adata: sc.AnnData,
    cluster_key: str = "clusters",
    **kwargs,
) -> None:
    """
    Neighborhood enrichment の結果をヒートマップとして可視化する。
    """
    sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, **kwargs)


def plot_co_occurrence(
    adata: sc.AnnData,
    cluster_key: str = "clusters",
    **kwargs,
) -> None:
    """
    Co-occurrence の結果をプロットする。
    """
    sq.pl.co_occurrence(adata, cluster_key=cluster_key, **kwargs)

