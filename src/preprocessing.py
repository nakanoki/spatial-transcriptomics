"""
前処理パイプライン（読み込み → QC → 正規化 → HVG）
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import scanpy as sc
import squidpy as sq


def read_visium_sample(sample_dir: Path) -> sc.AnnData:
    """
    Space Ranger 出力ディレクトリから Visium サンプルを読み込む。

    `*filtered_feature_bc_matrix.h5` と `spatial/` を含むディレクトリを想定する。
    10x の配布ファイルはサンプル名が接頭辞に付くため、glob で探索する。

    scanpy の `sc.read_visium` は squidpy へ移管され将来削除されるため、
    `sq.read.visium` を使う。
    """
    sample_dir = Path(sample_dir)
    h5_files = sorted(sample_dir.glob("*filtered_feature_bc_matrix.h5"))
    if not h5_files:
        raise FileNotFoundError(
            f"filtered_feature_bc_matrix.h5 が見つかりません: {sample_dir}"
        )
    adata = sq.read.visium(sample_dir, counts_file=h5_files[0].name)
    adata.var_names_make_unique()
    return adata


def load_processed_h5ad(path: Path) -> sc.AnnData:
    """処理済み .h5ad を読み込む"""
    return sc.read_h5ad(path)


def add_qc_metrics(adata: sc.AnnData, mt_prefix: str = "MT-") -> None:
    """
    ミトコンドリア遺伝子フラグと QC メトリクスを adata に追加する。
    """
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def filter_cells_and_genes(
    adata: sc.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    pct_counts_mt_max: float = 20.0,
) -> sc.AnnData:
    """
    セル・遺伝子のフィルタとミトコンドリア率でサブセット。コピーを返す。
    """
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < pct_counts_mt_max, :].copy()
    return adata


def normalize_and_hvg(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
) -> sc.AnnData:
    """
    正規化・log1p・高変動遺伝子選択。高変動遺伝子のみのコピーを返す。
    """
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    return adata



@dataclass(frozen=True)
class PreprocessParams:
    mt_prefix: str = "MT-"
    min_genes: int = 200
    min_cells: int = 3
    pct_counts_mt_max: float = 20.0
    target_sum: float = 1e4
    n_top_genes: int = 2000
    hvg_flavor: str = "seurat"
    normalization_method: Literal["log1p", "none"] = "log1p"


def preprocess_adata(
    adata: sc.AnnData,
    params: PreprocessParams = PreprocessParams(),
) -> sc.AnnData:
    """
    AnnData を前処理して返す。

    - QC 指標追加
    - セル/遺伝子フィルタ
    - 正規化 + log1p + HVG（HVG のみにサブセット）
    """
    adata = adata.copy()

    add_qc_metrics(adata, mt_prefix=params.mt_prefix)
    adata = filter_cells_and_genes(
        adata,
        min_genes=params.min_genes,
        min_cells=params.min_cells,
        pct_counts_mt_max=params.pct_counts_mt_max,
    )

    if params.normalization_method == "none":
        return adata

    adata = normalize_and_hvg(
        adata,
        target_sum=params.target_sum,
        n_top_genes=params.n_top_genes,
        flavor=params.hvg_flavor,
    )
    return adata

