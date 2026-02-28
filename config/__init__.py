"""
YAML 設定の読み込みユーティリティ。

Notebook からそのまま呼べるように、最小限の関数を提供します。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def load_config(project_root: Path, config_path: Path | None = None) -> dict[str, Any]:
    """
    `config/config.yaml` を dict として読み込む。
    """
    cfg_path = config_path or (project_root / "config" / "config.yaml")
    try:
        import yaml  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError("PyYAML が必要です（pip install pyyaml）") from e

    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    return data or {}


def get_paths(project_root: Path, cfg: dict[str, Any] | None = None) -> dict[str, Path]:
    """
    設定から主要パスを組み立てて返す。
    """
    cfg = cfg or load_config(project_root)
    paths_cfg = (cfg.get("paths") or {}) if isinstance(cfg, dict) else {}

    data_raw = Path(paths_cfg.get("data_raw_dir", "data/raw"))
    data_processed = Path(paths_cfg.get("data_processed_dir", "data/processed"))
    results = Path(paths_cfg.get("results_dir", "results"))
    notebooks = Path(paths_cfg.get("notebooks_dir", "notebooks"))
    scripts = Path(paths_cfg.get("scripts_dir", "scripts"))

    return {
        "project_root": project_root,
        "data_raw": (project_root / data_raw).resolve(),
        "data_processed": (project_root / data_processed).resolve(),
        "results": (project_root / results).resolve(),
        "notebooks": (project_root / notebooks).resolve(),
        "scripts": (project_root / scripts).resolve(),
    }


def get_qc_config(cfg: dict[str, Any] | None = None) -> dict[str, Any]:
    """
    Notebook 側が使いやすいキー名に寄せた QC 設定を返す。
    """
    cfg = cfg or {}
    qc_cfg = cfg.get("qc") or {}
    feat_cfg = cfg.get("feature_selection") or {}
    norm_cfg = cfg.get("normalization") or {}

    def _num(x, default):
        return default if x is None else x

    # config.yaml 側のキーと Notebook 既存コードのキーを折衷
    min_genes = _num(qc_cfg.get("min_genes_per_spot"), 0)
    min_cells = _num(qc_cfg.get("min_cells"), 3)  # 未定義なら一般的なデフォルト

    max_mito = qc_cfg.get("max_mito_percent")
    pct_counts_mt_max = 100.0 if max_mito is None else float(max_mito)

    return {
        "min_genes": int(min_genes),
        "min_cells": int(min_cells),
        "pct_counts_mt_max": float(pct_counts_mt_max),
        "n_top_genes": int(_num(feat_cfg.get("n_top_genes"), 2000)),
        "target_sum": float(_num(norm_cfg.get("target_sum"), 1e4)),
    }

