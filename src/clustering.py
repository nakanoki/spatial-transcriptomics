"""
次元削減・近傍グラフ・クラスタリング
"""

from __future__ import annotations

from dataclasses import dataclass

import scanpy as sc


def build_spatial_neighbors(
    adata: sc.AnnData, coord_type: str = "generic", delaunay: bool = True
) -> None:
    """
    Squidpy で空間隣接グラフを構築し adata に格納する。
    """
    try:
        import squidpy as sq  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError("squidpy が必要です（pip/conda で squidpy をインストールしてください）") from e

    sq.gr.spatial_neighbors(adata, coord_type=coord_type, delaunay=delaunay)


def run_moran_i(
    adata: sc.AnnData,
    genes: list[str] | None = None,
    use_highly_variable: bool = True,
    n_genes: int | None = 100,
) -> None:
    """
    モランI 空間自己相関を計算。結果は adata.uns["moranI"] に格納。
    """
    try:
        import squidpy as sq  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError("squidpy が必要です（pip/conda で squidpy をインストールしてください）") from e

    if genes is None and use_highly_variable and "highly_variable" in adata.var.columns:
        genes = adata.var_names[adata.var["highly_variable"]].tolist()
        if n_genes is not None:
            genes = genes[:n_genes]
    elif genes is None:
        genes = adata.var_names[: min(n_genes or 100, adata.n_vars)].tolist()

    sq.gr.spatial_autocorr(adata, mode="moran", genes=genes)


@dataclass(frozen=True)
class DimReduceParams:
    pca_n_components: int = 50
    neighbors_k: int = 15
    umap_min_dist: float = 0.5
    random_state: int = 0


@dataclass(frozen=True)
class ClusteringParams:
    method: str = "leiden"  # "leiden" or "louvain"
    resolution: float = 1.0
    key_added: str | None = None
    random_state: int = 0


def run_dimensionality_reduction(
    adata: sc.AnnData, params: DimReduceParams = DimReduceParams()
) -> sc.AnnData:
    """
    PCA → neighbors → UMAP を実行して返す（コピー）。
    """
    adata = adata.copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=params.pca_n_components, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=params.neighbors_k)
    sc.tl.umap(adata, min_dist=params.umap_min_dist, random_state=params.random_state)
    return adata


def run_clustering(
    adata: sc.AnnData, params: ClusteringParams = ClusteringParams()
) -> sc.AnnData:
    """
    Leiden/Louvain を実行して返す（コピー）。
    事前に neighbors が必要。
    """
    adata = adata.copy()

    method = params.method.lower()
    key_added = params.key_added or method

    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=params.resolution,
            key_added=key_added,
            random_state=params.random_state,
        )
        return adata

    if method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=params.resolution,
            key_added=key_added,
            random_state=params.random_state,
        )
        return adata

    raise ValueError(f"Unsupported clustering method: {params.method}")

