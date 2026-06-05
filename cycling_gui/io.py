from __future__ import annotations

import json
import re
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

import numpy as np

RESULTS_DIRNAME = "results"
ARCHIVE_NAME = "all_results.json"


def ensure_results_dir(base_dir: str | Path | None = None) -> Path:
    root = Path(base_dir) if base_dir is not None else Path.cwd()
    results_dir = root / RESULTS_DIRNAME
    results_dir.mkdir(parents=True, exist_ok=True)
    return results_dir


def safe_filename(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    return cleaned or "simulation"


def _jsonable(value: Any):
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.generic,)):
        return value.item()
    return value


def export_result_json(result, *, name: str, results_dir: str | Path | None = None) -> Path:
    results_dir = ensure_results_dir(results_dir)
    path = results_dir / f"{safe_filename(name)}.json"
    payload = result.to_dict() if hasattr(result, "to_dict") else _jsonable(result)
    path.write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")

    archive_path = results_dir / ARCHIVE_NAME
    archive = []
    if archive_path.exists():
        archive = json.loads(archive_path.read_text(encoding="utf-8"))
    archive.append(result.summary() if hasattr(result, "summary") else _jsonable(payload))
    archive_path.write_text(json.dumps(_jsonable(archive), indent=2), encoding="utf-8")
    return path

