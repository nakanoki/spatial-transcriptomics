"""
YAML 設定の読み込みユーティリティ。

Notebook からそのまま呼べるように、最小限の関数を提供する。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


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
