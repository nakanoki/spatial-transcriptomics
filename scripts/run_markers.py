"""
Issue #16: クラスタごとのマーカー遺伝子で空間ドメインの正体を同定する。

クラスタリング（scripts/run_clustering.py）の出力を入力とし、以下を出力する。

    results/markers/marker_genes.csv         クラスタごとの上位マーカー
    results/markers/known_marker_dotplot.png 既知マーカーの発現（dotplot）
    results/markers/top_marker_dotplot.png   各クラスタの上位マーカー
    results/markers/cluster_annotation.csv   クラスタごとの系統スコアと解釈
    results/markers/moran_cluster_map.csv    Moran's I 上位遺伝子とクラスタの対応

使い方:  uv run python scripts/run_markers.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 画面を持たない環境で実行するため
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import get_paths, load_config  # noqa: E402
from src.deg import DegParams, compute_markers, markers_to_dataframe  # noqa: E402
from src.preprocessing import load_processed_h5ad  # noqa: E402

CLUSTER_KEY = "clusters"

# 照合する既知マーカー。Visium のスポットは直径 55µm で複数細胞を含むため、
# 1スポットが単一の細胞種とは限らない。あくまで優勢な系統の手がかりとして使う。
KNOWN_MARKERS: dict[str, list[str]] = {
    "Epithelial/tumor": ["EPCAM", "KRT8", "KRT18", "ERBB2"],
    "Lymphocyte": ["PTPRC", "CD3D", "MS4A1"],
    "Plasma cell": ["IGKC", "IGHG1", "IGHG3", "MZB1", "JCHAIN"],
    "Macrophage": ["C1QA", "C1QB", "CD68"],
    "Stroma/fibroblast": ["COL1A1", "COL11A1", "DCN", "FAP"],
    "Endothelium": ["PECAM1", "VWF"],
    "Adipocyte": ["FABP4", "ADIPOQ", "PLIN1", "ADH1B"],
}


def present_markers(adata: sc.AnnData) -> dict[str, list[str]]:
    """データに存在する遺伝子だけに絞る。欠けているものは報告する。"""
    available = set(adata.raw.var_names)
    out, missing = {}, []
    for lineage, genes in KNOWN_MARKERS.items():
        found = [g for g in genes if g in available]
        missing += [g for g in genes if g not in available]
        if found:
            out[lineage] = found
    if missing:
        print(f"  データに存在しないマーカー: {', '.join(missing)}")
    return out


def score_lineages(adata: sc.AnnData, panels: dict[str, list[str]]) -> pd.DataFrame:
    """
    系統ごとのスコアをクラスタ平均で出す。

    `sc.tl.score_genes` は対照遺伝子群との比較でスコアを出すため、
    発現量の絶対値ではなく「その系統に偏っているか」を見られる。
    """
    for lineage, genes in panels.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_{lineage}", use_raw=True)
    cols = [f"score_{k}" for k in panels]
    scores = adata.obs.groupby(CLUSTER_KEY, observed=True)[cols].mean()
    scores.columns = list(panels)
    return scores


def zscore_by_lineage(scores: pd.DataFrame) -> pd.DataFrame:
    """
    系統ごとにクラスタ方向で z 化する。

    score_genes の値は「遺伝子セット平均 − 対照セット平均」で、パネル間で
    スケールが揃わない。乳腺組織では KRT8 / KRT18 / EPCAM がどこでも高いため、
    生の値で argmax を取ると大半のクラスタが上皮に寄ってしまう。
    系統ごとに z 化すれば「このクラスタはどの系統に相対的に偏っているか」を問える。
    """
    return (scores - scores.mean()) / scores.std(ddof=0)


def main() -> None:
    cfg = load_config(PROJECT_ROOT)
    paths = get_paths(PROJECT_ROOT, cfg=cfg)
    mk = cfg["markers"]

    out_dir = paths["results"] / "markers"
    out_dir.mkdir(parents=True, exist_ok=True)

    src_h5ad = paths["data_processed"] / f"{cfg['inputs']['sample_id']}_clustered.h5ad"
    if not src_h5ad.exists():
        raise FileNotFoundError(
            f"クラスタ済みデータが見つかりません: {src_h5ad}\n"
            "先に `uv run python scripts/run_clustering.py` を実行してください。"
        )
    adata = load_processed_h5ad(src_h5ad)
    if adata.raw is None:
        raise ValueError("adata.raw がありません。HVG subset 前の全遺伝子が必要です。")
    print(
        f"入力: {adata.n_obs:,} spots x {adata.n_vars:,} genes (HVG)"
        f" / adata.raw: {adata.raw.n_vars:,} genes"
    )

    # HVG 2,000 遺伝子では既知マーカーが落ちている可能性があるため raw を使う
    compute_markers(
        adata,
        DegParams(
            groupby=CLUSTER_KEY,
            method=str(mk["method"]),
            n_genes=int(mk["n_genes"]),
            use_raw=True,
        ),
    )
    markers = markers_to_dataframe(adata)
    markers.to_csv(out_dir / "marker_genes.csv", index=False)
    print(f"マーカー遺伝子: {len(markers):,} 行を保存")

    print("\n各クラスタの上位5マーカー:")
    for g, sub in markers.groupby("group", observed=True, sort=False):
        print(f"  cluster {g:>2}: {', '.join(sub['names'].head(5))}")

    # --- 既知マーカーとの照合 ---
    print("\n既知マーカーとの照合:")
    panels = present_markers(adata)
    scores = score_lineages(adata, panels)

    z = zscore_by_lineage(scores)

    sizes = adata.obs[CLUSTER_KEY].value_counts().sort_index()
    ann = z.copy()
    ann.insert(0, "spots", sizes)
    ann["assigned"] = z.idxmax(axis=1)
    # 最大と2番目の差。小さいほど判定が曖昧なので、断定を避ける目安にする
    top2 = z.apply(lambda r: r.nlargest(2).values, axis=1)
    ann["margin"] = [round(v[0] - v[1], 3) for v in top2]
    ann["confidence"] = [
        "high" if m >= 1.0 else ("medium" if m >= 0.4 else "low") for m in ann["margin"]
    ]
    scores.round(3).to_csv(out_dir / "lineage_scores_raw.csv", index_label="cluster")
    ann["mean_n_genes"] = (
        adata.obs.groupby(CLUSTER_KEY, observed=True)["n_genes_by_counts"].mean().round(0)
    )
    ann.round(3).to_csv(out_dir / "cluster_annotation.csv", index_label="cluster")
    print("（系統ごとにクラスタ方向で z 化した値）")
    print(ann.round(2).to_string())

    # --- Moran's I 上位遺伝子がどのクラスタのマーカーか ---
    if "moranI" in adata.uns:
        moran_top = adata.uns["moranI"].sort_values("I", ascending=False).head(15).index
        rows = []
        for gene in moran_top:
            if gene not in adata.raw.var_names:
                continue
            expr = pd.Series(
                adata.raw[:, gene].X.toarray().ravel(), index=adata.obs_names
            )
            per_cluster = expr.groupby(adata.obs[CLUSTER_KEY], observed=True).mean()
            rows.append(
                {
                    "gene": gene,
                    "moranI": round(float(adata.uns["moranI"].loc[gene, "I"]), 4),
                    "top_cluster": per_cluster.idxmax(),
                    "mean_expr": round(float(per_cluster.max()), 3),
                    "assigned_lineage": ann.loc[per_cluster.idxmax(), "assigned"],
                }
            )
        moran_map = pd.DataFrame(rows)
        moran_map.to_csv(out_dir / "moran_cluster_map.csv", index=False)
        print("\nMoran's I 上位遺伝子とクラスタの対応:")
        print(moran_map.to_string(index=False))

    # --- 図 ---
    sc.pl.dotplot(adata, panels, groupby=CLUSTER_KEY, use_raw=True, show=False)
    plt.savefig(out_dir / "known_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, show=False)
    plt.savefig(out_dir / "top_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    adata.write_h5ad(
        paths["data_processed"] / f"{cfg['inputs']['sample_id']}_annotated.h5ad"
    )
    print(f"\n保存: {out_dir}/")


if __name__ == "__main__":
    main()
