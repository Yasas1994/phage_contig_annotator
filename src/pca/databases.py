"""Database discovery and layout normalization for pca.

The expected database directory layout is hierarchical:

    databases/
    ├── phrogs/                       # database name
    │   ├── hmms/                     # HMM profiles (.hmm)
    │   │   ├── 9.hmm
    │   │   └── 99.hmm
    │   ├── metadata/                 # annotation metadata (.csv, .tsv)
    │   │   └── PHROG_annot_v4.csv
    │   ├── diamond/                  # DIAMOND database (.dmnd)
    │   └── mmseqs/                   # MMseqs database directory
    ├── apis/
    │   ├── hmms/
    │   │   └── dbAPIS.hmm
    │   └── metadata/
    │       └── seed_and_familyrep_all_infor.tsv
    └── db_chkpt

Each database is a top-level subdirectory of the database root. Inside a
database directory, specific subdirectories hold different database formats
or supporting files. The discovery code also tolerates legacy flat layouts
for backward compatibility.
"""
from __future__ import annotations

import logging
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

logger = logging.getLogger(__name__)

# Standard subdirectory names for different database formats.  Aliases are
# supported for discovery; normalization always produces the canonical names.
_HMM_SUBDIR = "hmms"
_HMM_ALIASES = {"hmms", "hmm"}
_METADATA_SUBDIR = "metadata"
_METADATA_ALIASES = {"metadata", "meta"}
_DIAMOND_SUBDIR = "diamond"
_MMSEQS_SUBDIR = "mmseqs"

# Recognised metadata file names (case-sensitive). The first match wins.
_METADATA_NAMES = ["PHROG_annot_v4.csv", "phrog_annot_v4.tsv"]


@dataclass
class DatabaseBundle:
    """A discovered database with optional format-specific subdirectories."""

    name: str
    root: Path
    hmms_path: Path | None = None
    metadata_path: Path | None = None
    diamond_path: Path | None = None
    mmseqs_path: Path | None = None

    @property
    def has_hmms(self) -> bool:
        return self.hmms_path is not None and self.hmms_path.exists()

    @property
    def has_metadata(self) -> bool:
        return self.metadata_path is not None and self.metadata_path.exists()

    @property
    def has_diamond(self) -> bool:
        return self.diamond_path is not None and self.diamond_path.exists()

    @property
    def has_mmseqs(self) -> bool:
        return self.mmseqs_path is not None and self.mmseqs_path.exists()


def discover_databases(db_dir: str | Path) -> dict[str, DatabaseBundle]:
    """Discover all databases under ``db_dir``.

    The function first tries the hierarchical layout (database subdirectories
    containing ``hmms/``, ``metadata/``, ``diamond/`` and ``mmseqs/``). If
    no database is found in the hierarchical layout, it falls back to the
    legacy flat layout where ``.hmm`` files and metadata files are discovered
    recursively.

    Returns
    -------
    dict[str, DatabaseBundle]
        Mapping from database name to its discovered bundle.
    """
    db_path = Path(db_dir)
    if not db_path.is_dir():
        return {}

    databases: dict[str, DatabaseBundle] = {}

    # Try hierarchical layout first.
    for subdir in sorted(p for p in db_path.iterdir() if p.is_dir()):
        bundle = _bundle_from_hierarchical_dir(subdir)
        if bundle is not None:
            databases[bundle.name] = bundle

    if databases:
        return databases

    # Fallback: legacy flat layout.
    logger.warning(
        "No hierarchical database layout found in %s; falling back to legacy discovery.",
        db_path,
    )
    return _discover_legacy_flat(db_path)


def _bundle_from_hierarchical_dir(root: Path) -> DatabaseBundle | None:
    """Return a bundle if ``root`` contains at least one standard subdirectory."""
    name = root.name

    # Find HMM subdirectory (canonical or alias).
    hmms_path = _find_subdir_by_name(root, _HMM_ALIASES)
    metadata = _find_metadata(root)
    diamond = _find_database_dir(root / _DIAMOND_SUBDIR)
    mmseqs = _find_database_dir(root / _MMSEQS_SUBDIR)

    valid = (
        hmms_path is not None
        or metadata is not None
        or diamond is not None
        or mmseqs is not None
    )
    if not valid:
        return None

    return DatabaseBundle(
        name=name,
        root=root,
        hmms_path=hmms_path,
        metadata_path=metadata,
        diamond_path=diamond,
        mmseqs_path=mmseqs,
    )


def _tsv_to_phrog_csv(tsv_path: Path, csv_path: Path) -> None:
    """Convert the legacy PHROG metadata TSV to the new CSV format.

    The new CSV keeps the same four core columns but uses CSV quoting, which
    is required for category names that contain commas (e.g. "DNA, RNA and
    nucleotide metabolism").
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.rename(
        columns={"phrog": "#phrog", "annot": "Annotation", "category": "Category"}
    )
    # The legacy TSV stores a numeric phrog id; the CSV stores the prefixed
    # identifier used by the search results.
    df["#phrog"] = df["#phrog"].apply(lambda x: f"phrog_{x}")
    df.to_csv(csv_path, index=False)


def _find_metadata(root: Path) -> Path | None:
    """Return the first recognised metadata file in a standard metadata dir.

    Checks both ``metadata/`` and the legacy ``meta/`` alias.
    """
    for subdir_name in _METADATA_ALIASES:
        metadata_dir = root / subdir_name
        if not metadata_dir.is_dir():
            continue
        for name in _METADATA_NAMES:
            candidate = metadata_dir / name
            if candidate.is_file():
                return candidate
    return None


def _find_database_dir(dbase_dir: Path) -> Path | None:
    """Return ``dbase_dir`` if it exists and is non-empty, else ``None``."""
    if not dbase_dir.is_dir():
        return None
    try:
        contents = list(dbase_dir.iterdir())
    except OSError:
        return None
    return dbase_dir if contents else None


def _find_subdir_by_name(root: Path, names: set[str]) -> Path | None:
    """Return the first existing subdirectory of ``root`` matching ``names``."""
    for name in names:
        candidate = root / name
        if candidate.is_dir():
            return candidate
    return None


def _discover_legacy_flat(db_path: Path) -> dict[str, DatabaseBundle]:
    """Discover databases from legacy flat or packaged layouts."""
    databases: dict[str, DatabaseBundle] = {}

    # Metadata discovery.
    metadata_path: Path | None = None
    for name in _METADATA_NAMES:
        candidates = list(db_path.rglob(name))
        if candidates:
            metadata_path = candidates[0]
            break

    # HMM discovery: group by parent directory.
    hmm_by_parent: dict[Path, list[Path]] = defaultdict(list)
    for hmm_file in db_path.rglob("*.hmm"):
        hmm_by_parent[hmm_file.parent].append(hmm_file)

    for parent, hmm_files in hmm_by_parent.items():
        parent_name = parent.name
        if parent == db_path:
            # The legacy flat layout (all HMMs loose at the root) is treated as
            # the default primary database for backward compatibility.
            name = "phrogs"
        elif "hmmerdb" in parent_name:
            name = parent_name.split("_")[0] if "_" in parent_name else "hmmerdb"
        else:
            name = parent.stem

        databases[name] = DatabaseBundle(
            name=name,
            root=parent,
            hmms_path=parent,
            metadata_path=metadata_path,
        )

    return databases


def normalize_database_dir(db_dir: str | Path) -> dict[str, DatabaseBundle]:
    """Organize files under ``db_dir`` into the hierarchical layout.

    Files are grouped by database name (top-level subdirectory) and then moved
    into the appropriate format subdirectories. Existing hierarchical layouts
    are preserved. The function is idempotent: running it twice on the same
    directory does not change the result.
    """
    db_path = Path(db_dir)
    db_path.mkdir(parents=True, exist_ok=True)

    # Phase 1: ensure already-hierarchical databases are represented.
    databases = discover_databases(db_path)

    # Phase 1b: rename legacy alias subdirectories to canonical names.
    _consolidate_alias_dirs(db_path, databases)
    databases = discover_databases(db_path)

    # Phase 2: move loose HMM files into database directories.
    loose_hmms = _find_loose_hmms(db_path, databases)
    _group_files_by_name(loose_hmms, _infer_database_name, _HMM_SUBDIR, db_path)

    # Phase 3: move loose metadata files into their database metadata dirs.
    loose_metadata = _find_loose_metadata(db_path, databases)
    if loose_metadata:
        metadata_by_db: dict[str, list[Path]] = defaultdict(list)
        for meta_file in loose_metadata:
            metadata_by_db[_infer_database_name(meta_file, db_path)].append(meta_file)
        for target_db, meta_files in metadata_by_db.items():
            target_dir = db_path / target_db / _METADATA_SUBDIR
            target_dir.mkdir(parents=True, exist_ok=True)
            for meta_file in meta_files:
                # Convert legacy PHROG TSV to CSV on the fly.
                if meta_file.name == "phrog_annot_v4.tsv":
                    dest = target_dir / "PHROG_annot_v4.csv"
                    _tsv_to_phrog_csv(meta_file, dest)
                    if meta_file.parent == db_path:
                        meta_file.unlink()
                    continue

                dest = target_dir / meta_file.name
                if dest != meta_file:
                    if dest.exists():
                        shutil.copy2(meta_file, dest)
                        meta_file.unlink()
                    else:
                        meta_file.rename(dest)

    # Phase 4: move loose DIAMOND files.
    loose_dmnd = _find_loose_files(db_path, "*.dmnd", databases)
    _group_files_by_name(loose_dmnd, _infer_database_name, _DIAMOND_SUBDIR, db_path)

    # Phase 5: clean up empty directories left behind.
    _remove_empty_dirs(db_path)

    return discover_databases(db_path)


def _consolidate_alias_dirs(
    db_path: Path, databases: dict[str, DatabaseBundle]
) -> None:
    """Move files from ``hmm/`` -> ``hmms/`` and ``meta/`` -> ``metadata/``."""
    for subdir in sorted(p for p in db_path.iterdir() if p.is_dir()):
        if subdir.name in databases and databases[subdir.name].hmms_path is not None:
            # Already handled by the canonical bundle; nothing to do.
            pass

        hmm_alias = _find_subdir_by_name(subdir, _HMM_ALIASES - {_HMM_SUBDIR})
        canonical_hmm = subdir / _HMM_SUBDIR
        if hmm_alias is not None and not canonical_hmm.exists():
            hmm_alias.rename(canonical_hmm)

        meta_alias = _find_subdir_by_name(subdir, _METADATA_ALIASES - {_METADATA_SUBDIR})
        canonical_meta = subdir / _METADATA_SUBDIR
        if meta_alias is not None and not canonical_meta.exists():
            meta_alias.rename(canonical_meta)


def _find_loose_hmms(
    db_path: Path, databases: dict[str, DatabaseBundle]
) -> list[Path]:
    """Return HMM profiles and index files not inside a standard HMM directory."""
    protected = set()
    for bundle in databases.values():
        if bundle.hmms_path is not None and bundle.hmms_path.name in _HMM_ALIASES:
            protected.add(bundle.hmms_path.resolve())

    hmm_patterns = ("*.hmm", "*.h3f", "*.h3i", "*.h3m", "*.h3p")
    loose: list[Path] = []
    for pattern in hmm_patterns:
        for hmm_file in db_path.rglob(pattern):
            if hmm_file.parent.resolve() in protected:
                continue
            if any(part in _HMM_ALIASES for part in hmm_file.parts):
                continue
            loose.append(hmm_file)
    return loose


def _find_loose_metadata(
    db_path: Path, databases: dict[str, DatabaseBundle]
) -> list[Path]:
    """Return metadata files not already inside a standard metadata directory."""
    protected: set[Path] = set()
    for bundle in databases.values():
        if bundle.metadata_path is not None and bundle.metadata_path.parent.name in _METADATA_ALIASES:
            protected.add(bundle.metadata_path.parent.resolve())

    loose: list[Path] = []
    for name in _METADATA_NAMES:
        for meta_file in db_path.rglob(name):
            if meta_file.parent.resolve() in protected:
                continue
            if any(part in _METADATA_ALIASES for part in meta_file.parts):
                continue
            loose.append(meta_file)
    return loose


def _find_loose_files(
    db_path: Path, pattern: str, databases: dict[str, DatabaseBundle]
) -> list[Path]:
    """Return files matching ``pattern`` not inside protected standard dirs."""
    protected: set[Path] = set()
    for bundle in databases.values():
        if bundle.diamond_path is not None and bundle.diamond_path.name == _DIAMOND_SUBDIR:
            protected.add(bundle.diamond_path.resolve())
        if bundle.mmseqs_path is not None and bundle.mmseqs_path.name == _MMSEQS_SUBDIR:
            protected.add(bundle.mmseqs_path.resolve())

    loose: list[Path] = []
    for file in db_path.rglob(pattern):
        if file.parent.resolve() in protected:
            continue
        loose.append(file)
    return loose


def _group_files_by_name(
    files: Iterable[Path],
    name_func: Callable[[Path, Path], str],
    subdir_name: str,
    db_path: Path,
) -> None:
    """Move ``files`` into ``<db_path>/<name>/<subdir_name>/``.

    ``name_func`` receives the file path and the database root path and
    returns the database name.
    """
    by_name: dict[str, list[Path]] = defaultdict(list)
    for file in files:
        by_name[name_func(file, db_path)].append(file)

    for name, file_list in by_name.items():
        target_dir = db_path / name / subdir_name
        target_dir.mkdir(parents=True, exist_ok=True)
        for file in file_list:
            dest = target_dir / file.name
            if dest != file:
                if dest.exists():
                    shutil.copy2(file, dest)
                    file.unlink()
                else:
                    file.rename(dest)


def _infer_database_name(file: Path, db_path: Path) -> str:
    """Infer a database name from a file path.

    If the file sits directly under ``db_path``, default to ``phrogs`` (or the
    only existing database subdirectory). Otherwise use the immediate
    subdirectory name under ``db_path``.
    """
    try:
        rel = file.relative_to(db_path)
    except ValueError:
        return "phrogs"
    if len(rel.parts) == 1:
        # File is directly under db_path; prefer phrogs if no other hint.
        return "phrogs"
    return rel.parts[0]


def _remove_empty_dirs(root: Path) -> None:
    """Remove empty directories under ``root`` except the root itself."""
    for subdir in sorted(root.rglob("*"), reverse=True):
        if subdir.is_dir():
            try:
                subdir.rmdir()
            except OSError:
                pass


# Need to import Callable here to avoid forward-reference issues.
from typing import Callable  # noqa: E402
