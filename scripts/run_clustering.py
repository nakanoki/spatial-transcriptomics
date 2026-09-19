"""
Issue #4: クラスタリングと空間ドメインの可視化。

前処理（scripts/run_preprocess.py）の出力を入力とし、以下を出力する。

    results/clustering/umap_clusters.png       クラスタで着色した UMAP
    results/clustering/spatial_clusters.png    組織像上の空間ドメイン
    results/clustering/moran_top_genes.csv     空間自己相関の上位遺伝子
    results/clustering/cluster_sizes.csv       各クラスタのスポット数
    data/processed/<sample_id>_clustered.h5ad  クラスタ付きデータ

使い方:  uv run python scripts/run_clustering.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 画面を持たない環境で実行するため

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import get_paths, load_config  # noqa: E402
from src.clustering import (  # noqa: E402
    ClusteringParams,
    DimReduceParams,
    build_spatial_neighbors,
    run_clustering,
    run_dimensionality_reduction,
    run_moran_i,
)
from src.preprocessing import load_processed_h5ad  # noqa: E402
from src.visualization import plot_spatial, plot_umap  # noqa: E402

CLUSTER_KEY = "clusters"
# resolution の妥当性を見るために比較する候補
SWEEP = (0.5, 0.8, 1.0, 1.5)


def sweep_resolution(adata, method: str, seed: int) -> None:
    """
    resolution を変えたときのクラスタ数を並べて出す。

    値を決めるための材料であり、結果は保存しない。近傍グラフは共通なので
    Leiden を回し直すだけで済む。
    """
    print("\nresolution ごとのクラスタ数:")
    for res in SWEEP:
        tmp = run_clustering(
            adata,
            ClusteringParams(
                method=method, resolution=res, key_added=CLUSTER_KEY, random_state=seed
            ),
        )
        sizes = tmp.obs[CLUSTER_KEY].value_counts()
        tiny = int((sizes < 50).sum())
        print(
            f"  resolution={res:<4} クラスタ数 {sizes.size:>3}"
            f"  最小 {sizes.min():>4} spots"
            f"  50 spots 未満のクラスタ {tiny}"
        )


def main() -> None:
    cfg = load_config(PROJECT_ROOT)
    paths = get_paths(PROJECT_ROOT, cfg=cfg)

    seed = int((cfg.get("project") or {}).get("random_seed") or 0)
    dr = cfg["dimension_reduction"]
    cl = cfg["clustering"]

    out_dir = paths["results"] / "clustering"
    out_dir.mkdir(parents=True, exist_ok=True)

    src_h5ad = paths["data_processed"] / f"{cfg['inputs']['sample_id']}_preprocessed.h5ad"
    if not src_h5ad.exists():
        raise FileNotFoundError(
            f"前処理済みデータが見つかりません: {src_h5ad}\n"
            "先に `uv run python scripts/run_preprocess.py` を実行してください。"
        )
    adata = load_processed_h5ad(src_h5ad)
    print(f"入力（前処理済み）: {adata.n_obs:,} spots x {adata.n_vars:,} genes")

    adata = run_dimensionality_reduction(
        adata,
        DimReduceParams(
            pca_n_components=int(dr["pca_n_components"]),
            neighbors_k=int(dr["neighbors_k"]),
            umap_min_dist=float(dr["umap_min_dist"]),
            random_state=seed,
        ),
    )

    sweep_resolution(adata, str(cl["method"]), seed)

    resolution = float(cl["resolution"])
    adata = run_clustering(
        adata,
        ClusteringParams(
            method=str(cl["method"]),
            resolution=resolution,
            key_added=CLUSTER_KEY,
            random_state=seed,
        ),
    )
    sizes = adata.obs[CLUSTER_KEY].value_counts().sort_index()
    print(f"\n採用 resolution={resolution} → {sizes.size} クラスタ")
    print(sizes.to_string())
    sizes.rename("spots").to_csv(out_dir / "cluster_sizes.csv", index_label="cluster")

    # 空間自己相関。空間的に構造を持つ遺伝子を確認する
    build_spatial_neighbors(adata, coord_type="generic", delaunay=True)
    run_moran_i(adata, n_genes=100)
    moran = adata.uns["moranI"].sort_values("I", ascending=False)
    moran.head(30).to_csv(out_dir / "moran_top_genes.csv", index_label="gene")
    print("\nMoran's I 上位 10:")
    print(moran.head(10)[["I", "pval_norm_fdr_bh"]].round(4).to_string())

    label = f"{cfg['inputs']['sample_id']}  |  {cl['method']} resolution={resolution}  |  {sizes.size} clusters"
    plot_umap(
        adata, color=CLUSTER_KEY,
        save=out_dir / "umap_clusters.png",
        title=f"UMAP colored by cluster  |  {label}",
    )
    plot_spatial(
        adata, color=CLUSTER_KEY,
        save=out_dir / "spatial_clusters.png",
        title=f"Spatial domains  |  {label}",
    )

    out_h5ad = paths["data_processed"] / f"{cfg['inputs']['sample_id']}_clustered.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"\n保存:\n  {out_h5ad}\n  {out_dir}/")


if __name__ == "__main__":
    main()
