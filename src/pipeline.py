"""
設定ファイル（config/config.yaml）に基づく簡易パイプライン。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import scanpy as sc

from config import get_paths, load_config
from .clustering import (
    ClusteringParams,
    DimReduceParams,
    build_spatial_neighbors,
    run_clustering,
    run_dimensionality_reduction,
    run_moran_i,
)
from .preprocessing import PreprocessParams, preprocess_adata, read_visium_sample


@dataclass(frozen=True)
class PipelineOutputs:
    adata: sc.AnnData
    output_dir: Path


def _pick_visium_input_dir(project_root: Path, cfg: dict[str, Any], paths: dict[str, Path]) -> Path:
    inputs = cfg.get("inputs") or {}
    spaceranger_out_dir = (inputs.get("spaceranger_out_dir") or "").strip()
    if spaceranger_out_dir:
        return (project_root / spaceranger_out_dir).resolve()

    sample_id = (inputs.get("sample_id") or "").strip()
    if sample_id:
        return (paths["data_raw"] / sample_id).resolve()

    samples = inputs.get("samples") or []
    if isinstance(samples, list) and samples:
        first = samples[0] or {}
        p = (first.get("spaceranger_out_dir") or "").strip()
        if p:
            return (project_root / p).resolve()
        sid = (first.get("sample_id") or "").strip()
        if sid:
            return (paths["data_raw"] / sid).resolve()

    raise ValueError(
        "入力データの場所が未設定です。config.yaml の inputs.spaceranger_out_dir "
        "または inputs.sample_id を設定してください。"
    )


def run_pipeline(project_root: Path, config_path: Path | None = None) -> PipelineOutputs:
    """
    Visium を想定したパイプラインを実行して、AnnData と出力先ディレクトリを返す。
    """
    cfg = load_config(project_root, config_path=config_path)
    paths = get_paths(project_root, cfg=cfg)

    # 入力決定
    visium_dir = _pick_visium_input_dir(project_root, cfg, paths)

    # 乱数シード
    seed = int(((cfg.get("project") or {}).get("random_seed") or 0))

    # 前処理
    qc_cfg = cfg.get("qc") or {}
    norm_cfg = cfg.get("normalization") or {}
    feat_cfg = cfg.get("feature_selection") or {}

    preprocess_params = PreprocessParams(
        mt_prefix="MT-",
        min_genes=int(qc_cfg.get("min_genes_per_spot") or 0),
        min_cells=int(qc_cfg.get("min_cells") or 3),
        pct_counts_mt_max=float(qc_cfg.get("max_mito_percent") or 100.0),
        target_sum=float(norm_cfg.get("target_sum") or 1e4),
        n_top_genes=int(feat_cfg.get("n_top_genes") or 2000),
        normalization_method=("none" if (norm_cfg.get("method") == "none") else "log1p"),
    )

    adata = read_visium_sample(visium_dir)
    adata = preprocess_adata(adata, preprocess_params)

    # 次元削減・クラスタリング
    dr_cfg = cfg.get("dimension_reduction") or {}
    cl_cfg = cfg.get("clustering") or {}

    adata = run_dimensionality_reduction(
        adata,
        DimReduceParams(
            pca_n_components=int(dr_cfg.get("pca_n_components") or 50),
            neighbors_k=int(dr_cfg.get("neighbors_k") or 15),
            umap_min_dist=float(dr_cfg.get("umap_min_dist") or 0.5),
            random_state=seed,
        ),
    )

    cluster_key = "clusters"
    adata = run_clustering(
        adata,
        ClusteringParams(
            method=str(cl_cfg.get("method") or "leiden"),
            resolution=float(cl_cfg.get("resolution") or 1.0),
            key_added=cluster_key,
            random_state=seed,
        ),
    )

    # 空間解析（squidpy は任意）
    try:
        build_spatial_neighbors(adata, coord_type="generic", delaunay=True)
        run_moran_i(adata, n_genes=100)
    except Exception:
        # squidpy 未導入などは許容（パイプライン全体は止めない）
        pass

    # 出力
    export_cfg = cfg.get("export") or {}
    run_name = str(export_cfg.get("run_name") or "run")
    out_dir = (paths["results"] / run_name).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if bool(export_cfg.get("save_anndata", True)):
        adata.write_h5ad(out_dir / "adata.h5ad")

    if bool(export_cfg.get("save_plots", True)):
        import matplotlib.pyplot as plt

        dpi = int(export_cfg.get("dpi") or 150)

        # UMAP
        sc.pl.umap(adata, color=[cluster_key], show=False)
        plt.savefig(out_dir / "umap_clusters.png", dpi=dpi, bbox_inches="tight")
        plt.close()

        # Spatial（Visium の場合）
        if "spatial" in adata.obsm:
            try:
                sc.pl.spatial(adata, color=[cluster_key], alpha_img=0.8, show=False)
                plt.savefig(out_dir / "spatial_clusters.png", dpi=dpi, bbox_inches="tight")
                plt.close()
            except Exception:
                pass

    return PipelineOutputs(adata=adata, output_dir=out_dir)

