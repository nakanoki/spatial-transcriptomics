"""
空間トランスクリプトームにおける差次的発現（DEG / マーカー遺伝子）解析ヘルパー。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd
import scanpy as sc


@dataclass(frozen=True)
class DegParams:
    groupby: str
    method: str = "wilcoxon"  # 't-test', 't-test_overestim_var', 'wilcoxon', 'logreg' など
    n_genes: int | None = None
    key_added: str = "rank_genes_groups"
    use_raw: bool | None = None
    layer: str | None = None


def compute_markers(
    adata: sc.AnnData,
    params: DegParams,
    **kwargs: Any,
) -> None:
    """
    `sc.tl.rank_genes_groups` をラップしてマーカー遺伝子を計算する。

    結果は `adata.uns[params.key_added]` に格納される。
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=params.groupby,
        method=params.method,
        n_genes=params.n_genes,
        key_added=params.key_added,
        use_raw=params.use_raw,
        layer=params.layer,
        **kwargs,
    )


def markers_to_dataframe(
    adata: sc.AnnData,
    group: str | None = None,
    key: str = "rank_genes_groups",
) -> pd.DataFrame:
    """
    `sc.get.rank_genes_groups_df` を使って、マーカー遺伝子結果を DataFrame に変換する。

    Parameters
    ----------
    adata:
        `compute_markers` を実行済みの AnnData
    group:
        特定のクラスター/ラベル名。None の場合は全グループ分を結合した DataFrame を返す。
    key:
        `adata.uns` に格納された結果のキー
    """
    df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
    # 一般的によく使う列名だけに絞る例（必要に応じてカスタマイズ）
    keep_cols = [c for c in df.columns if c in {"group", "names", "scores", "logfoldchanges", "pvals_adj"}]
    if keep_cols:
        df = df[keep_cols]
    return df

