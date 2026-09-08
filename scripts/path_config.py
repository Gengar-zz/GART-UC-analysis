from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RepositoryPaths:
    repo_root: Path
    external_data_root: Path
    work_root: Path
    data_public_root: Path
    results_root: Path


def get_paths() -> RepositoryPaths:
    repo_root = Path(__file__).resolve().parents[1]
    override = os.environ.get("GART_UC_DATA_ROOT", "").strip()
    external = Path(override).expanduser().resolve() if override else repo_root / "data_external"
    return RepositoryPaths(repo_root, external, repo_root / "work", repo_root / "data_public", repo_root / "results_reproduced")


def require_file(path: Path, label: str = "required input") -> Path:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Required external input not found: {path} ({label}). See data_external/README.md or set GART_UC_DATA_ROOT.")
    return path
