"""Command-line interface for the phage contig annotator pipeline."""

from __future__ import annotations

import configparser
import os
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Sequence

import click
import requests
import tqdm

from pca import annotations, validation
from pca.logutils import get_logger


def _package_dir() -> Path:
    """Return the directory containing this package."""
    return Path(__file__).parent.resolve()


_DEFAULT_DB_DIR: Path = _package_dir().parents[1] / "databases"


def _load_config() -> configparser.ConfigParser:
    """Load the bundled ``config.ini``."""
    config_path = _package_dir() / "data" / "config.ini"
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def _write_config(config: configparser.ConfigParser) -> None:
    """Persist updates to the bundled ``config.ini``."""
    config_path = _package_dir() / "data" / "config.ini"
    with open(config_path, "w") as configfile:
        config.write(configfile)


def _tsv_to_phrog_csv(tsv_path: Path, csv_path: Path) -> None:
    """Convert the legacy PHROG metadata TSV to the new CSV format.

    The new CSV keeps the same four core columns but uses CSV quoting, which
    is required for category names that contain commas (e.g. "DNA, RNA and
    nucleotide metabolism").
    """
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t")
    df = df.rename(
        columns={"phrog": "#phrog", "annot": "Annotation", "category": "Category"}
    )
    # The legacy TSV stores a numeric phrog id; the CSV stores the prefixed
    # identifier used by the search results.
    df["#phrog"] = df["#phrog"].apply(lambda x: f"phrog_{x}")
    df.to_csv(csv_path, index=False)


def _discover_databases(db_dir: str) -> tuple[dict[str, str], str | None]:
    """Return discovered HMM databases and the PHROG metadata path.

    Supports both the expected packaged layout (``hmmerdb*`` directories and a
    ``meta`` directory containing ``PHROG_annot_v4.csv``) and flat layouts
    produced by older extraction code.
    """
    db_path = Path(db_dir)
    hmmerdb: dict[str, str] = {}
    meta_path: str | None = None

    # Prefer the extended PHROG metadata CSV; fall back to the legacy TSV.
    meta_candidates = list(db_path.rglob("PHROG_annot_v4.csv"))
    if not meta_candidates:
        meta_candidates = list(db_path.rglob("phrog_annot_v4.tsv"))
    if meta_candidates:
        meta_path = str(meta_candidates[0])

    # Locate HMM files and group them by their parent directory.
    hmm_paths = list(db_path.rglob("*.hmm"))
    if not hmm_paths:
        return hmmerdb, meta_path

    for hmm_file in hmm_paths:
        parent = hmm_file.parent
        parent_name = parent.name
        # Prefer directories that look like the packaged hmmerdb splits.
        if "hmmerdb" in parent_name:
            db_name = parent_name.split("_")[0]
        elif parent == db_path:
            db_name = "hmmerdb"
        else:
            db_name = parent.stem
        hmmerdb[db_name] = str(parent)

    return hmmerdb, meta_path


def _snakefile_path() -> Path:
    """Return the path to the bundled Snakemake workflow."""
    return _package_dir() / "workflow" / "Snakefile"


def download_dbs(path: str) -> bool:
    """Download and extract the PHROG annotation database."""
    checkpoint = Path(path) / "db_chkpt"
    if checkpoint.exists():
        click.echo(f"Skipping database download as checkpoint found at {path}")
        return True

    url = (
        "https://nextcloud.uni-greifswald.de/index.php/s/"
        "ft8FAoQXscoj9eo/download/database.tar.gz"
    )
    os.makedirs(path, exist_ok=True)

    response = requests.get(url, stream=True, timeout=300)
    if response.status_code != 200:
        click.echo(
            f"Failed to download the database from {url}. "
            f"Status code: {response.status_code}",
            err=True,
        )
        return False

    tar_file_path = Path(path) / "database.tar.gz"
    total_size = int(response.headers.get("content-length", 0))
    pbar = tqdm.tqdm(total=total_size / (1024 * 1024), unit="MB")
    with open(tar_file_path, "wb") as file:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                pbar.update(len(chunk) / (1024 * 1024))
                file.write(chunk)
    pbar.close()

    base_path = Path(path).resolve()
    with tarfile.open(tar_file_path, "r:gz") as tar:
        for member in tar.getmembers():
            # Strip any leading path components and reject unsafe members.
            member_path = Path(member.name)
            if not member.isfile():
                continue
            if any(part == ".." for part in member_path.parts):
                continue
            if member_path.is_absolute():
                continue
            # Resolve relative to the target directory to prevent traversal.
            target = (base_path / member_path).resolve()
            if not str(target).startswith(str(base_path)):
                continue
            member.name = str(member_path)
        tar.extractall(path=path)

    _normalize_database_dir(path)

    tar_file_path.unlink()
    click.echo(f"Database downloaded and extracted to {path}")
    checkpoint.touch()
    return True


def _normalize_database_dir(path: str) -> None:
    """Organize extracted database files into a predictable layout.

    Moves loose ``.hmm`` files into a ``hmmerdb/`` directory and the PHROG
    metadata CSV into a ``meta/`` directory when they are not already there.
    """
    db_path = Path(path)
    hmmerdb_dir = db_path / "hmmerdb"
    meta_dir = db_path / "meta"

    # Gather loose .hmm files (those not already under a hmmerdb* directory).
    loose_hmms = [
        p for p in db_path.rglob("*.hmm")
        if "hmmerdb" not in p.parent.name
    ]
    if loose_hmms:
        hmmerdb_dir.mkdir(exist_ok=True)
        for hmm_file in loose_hmms:
            dest = hmmerdb_dir / hmm_file.name
            if dest != hmm_file:
                hmm_file.rename(dest)

    # Gather loose metadata CSVs.
    csv_files = [p for p in db_path.rglob("PHROG_annot_v4.csv") if p.parent != meta_dir]
    if csv_files:
        meta_dir.mkdir(exist_ok=True)
        for meta_file in csv_files:
            dest = meta_dir / meta_file.name
            if dest != meta_file:
                meta_file.rename(dest)

    # Convert any loose legacy TSVs to the new CSV format when no CSV is present.
    tsv_files = [
        p for p in db_path.rglob("phrog_annot_v4.tsv")
        if p.parent != meta_dir and not (meta_dir / "PHROG_annot_v4.csv").exists()
    ]
    if tsv_files:
        meta_dir.mkdir(exist_ok=True)
        for tsv_file in tsv_files:
            csv_dest = meta_dir / "PHROG_annot_v4.csv"
            _tsv_to_phrog_csv(tsv_file, csv_dest)
            if tsv_file.parent == db_path:
                tsv_file.unlink()

    # Remove empty directories left behind after moving files.
    for subdir in sorted(db_path.iterdir(), reverse=True):
        if subdir.is_dir() and subdir not in (hmmerdb_dir, meta_dir):
            try:
                subdir.rmdir()
            except OSError:
                pass


@click.group(invoke_without_command=True)
@click.option("--quiet", is_flag=True, help="Suppress informational logging.")
@click.pass_context
def main(ctx: click.Context, quiet: bool) -> None:
    """Annotate genes in phage contigs with PHROG/VOG HMMs."""
    ctx.ensure_object(dict)
    ctx.obj["quiet"] = quiet
    get_logger(quiet=quiet)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)


@main.command()
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Path to input FASTA file with contigs or proteins.",
)
@click.option(
    "-o",
    "--output",
    "output_dir",
    required=True,
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    help="Path to output directory.",
)
@click.option(
    "--type",
    "input_type",
    default="contigs",
    type=click.Choice(["contigs", "proteins"], case_sensitive=False),
    help="Input type (default: contigs).",
)
@click.option(
    "--skip-trna",
    is_flag=True,
    help="Skip tRNAscan-SE gene prediction.",
)
@click.option(
    "-db",
    "--db",
    "db_dir",
    default=None,
    type=click.Path(exists=True, file_okay=False, readable=True, path_type=Path),
    help="Use a custom HMM database directory.",
)
@click.option(
    "--cpus",
    default=8,
    show_default=True,
    type=int,
    help="Number of CPUs to use.",
)
@click.option(
    "--plot-format",
    "plot_formats",
    multiple=True,
    default=("pdf",),
    show_default=True,
    type=click.Choice(["pdf", "png", "html"], case_sensitive=False),
    help="Plot output format. Specify multiple times for multiple formats.",
)
@click.option(
    "--translation-table",
    "-t",
    default=None,
    type=int,
    help="NCBI translation table to force for gene calling (e.g. 11 or 4). "
         "By default pyrodigal-gv selects the table automatically in metagenomic mode.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show Snakemake execution plan without running.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Force re-execution of all workflow steps.",
)
@click.option(
    "--theme",
    default="light",
    show_default=True,
    type=click.Choice(["light", "dark"], case_sensitive=False),
    help="Color theme for generated plots.",
)
@click.pass_context
def run(
    ctx: click.Context,
    input_path: Path,
    output_dir: Path,
    input_type: str,
    skip_trna: bool,
    db_dir: Path | None,
    cpus: int,
    plot_formats: tuple[str, ...],
    translation_table: int | None,
    dry_run: bool,
    force: bool,
    theme: str,
) -> None:
    """Run the full annotation pipeline."""
    logger = get_logger(quiet=ctx.obj["quiet"])

    config = _load_config()
    db_dir = db_dir or Path(config["databases"]["dbroot"] or _DEFAULT_DB_DIR)
    if not db_dir.is_dir():
        raise click.UsageError(f"Database directory does not exist: {db_dir}")

    hmmerdb, meta_path = _discover_databases(str(db_dir))
    if not hmmerdb:
        raise click.UsageError(
            f"No HMM databases found in {db_dir}. "
            "Run 'phage_contig_annotator download-db' first."
        )
    if meta_path is None:
        raise click.UsageError(
            f"PHROG metadata not found in {db_dir}. "
            "Run 'phage_contig_annotator download-db' first."
        )

    if not skip_trna:
        validation.check_executables(["tRNAscan-SE"])

    output_dir.mkdir(parents=True, exist_ok=True)

    workflow_config = {
        "input": str(input_path.resolve()),
        "input_type": input_type,
        "output_dir": str(output_dir.resolve()),
        "db_dir": str(db_dir.resolve()),
        "meta_path": str(Path(meta_path).resolve()),
        "run_trna": not skip_trna,
        "plot_formats": list(plot_formats),
        "translation_table": translation_table,
        "theme": theme,
    }

    config_path = output_dir / "config.yaml"
    _write_yaml_config(config_path, workflow_config)
    logger.info("wrote workflow config to %s", config_path)

    snakefile = _snakefile_path()
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--directory", str(output_dir),
        "--cores", str(cpus),
        "--configfile", str(config_path),
    ]
    if dry_run:
        cmd.append("--dry-run")
    if force:
        cmd.append("--forceall")
        cmd.append("--rerun-incomplete")

    env = os.environ.copy()
    src_dir = str(_package_dir().parent)
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{src_dir}{os.pathsep}{existing_pythonpath}" if existing_pythonpath else src_dir

    cmd.extend(["--envvars", "PYTHONPATH"])

    if force:
        unlock_cmd = cmd + ["--unlock"]
        logger.info("unlocking snakemake directory: %s", " ".join(unlock_cmd))
        subprocess.run(unlock_cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    logger.info("running snakemake: %s", " ".join(cmd))
    result = subprocess.run(cmd, env=env)
    sys.exit(result.returncode)


@main.command()
@click.option(
    "-p",
    "--path",
    "-db",
    "--db",
    "path",
    default=None,
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    help="Path to the database directory (alias: --db).",
)
@click.pass_context
def download_db(ctx: click.Context, path: Path | None) -> None:
    """Download the annotation database and update local config."""
    path = path or _DEFAULT_DB_DIR
    path.mkdir(parents=True, exist_ok=True)

    if not download_dbs(path=str(path)):
        raise click.ClickException("Database download failed.")

    config = _load_config()
    config.set("databases", "dbroot", str(path))
    config.set("databases", "hmmdb", str(path / "hmmdb"))
    config.set("databases", "meta", str(path / "meta"))
    _write_config(config)
    click.echo(f"Updated config.ini with database path: {path}")


@main.group()
@click.pass_context
def utils(ctx: click.Context) -> None:
    """Utility commands for the annotator."""


@utils.command("report")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Path to a GenBank (.gb, .gbk) or GFF (.gff, .gff3) file.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(file_okay=True, dir_okay=True, writable=True, path_type=Path),
    help="Output .html file (single contig) or output directory (multiple contigs).",
)
@click.option(
    "-f",
    "--fasta",
    "fasta_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Optional nucleotide FASTA file for GFF inputs without embedded sequences.",
)
@click.option(
    "--meta",
    "meta_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Optional PHROG metadata CSV for category colors.",
)
@click.option(
    "--theme",
    default="light",
    show_default=True,
    type=click.Choice(["light", "dark"], case_sensitive=False),
    help="Color theme for the HTML report.",
)
@click.pass_context
def utils_report(
    ctx: click.Context,
    input_path: Path,
    output_path: Path,
    fasta_path: Path | None,
    meta_path: Path | None,
    theme: str,
) -> None:
    """Convert a GenBank or GFF file to an interactive HTML report."""
    logger = get_logger(quiet=ctx.obj["quiet"])
    try:
        out_paths = annotations.convert_to_html(
            input_path,
            output_path,
            fasta_path=fasta_path,
            meta_path=meta_path,
            theme=theme,
        )
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    for out_path in out_paths:
        click.echo(f"Wrote {out_path}")
    logger.info("generated %d HTML report(s)", len(out_paths))


def _write_yaml_config(path: Path, config: dict[str, object]) -> None:
    """Write a workflow config YAML file."""
    import yaml

    with open(path, "w") as fh:
        yaml.safe_dump(config, fh, default_flow_style=False)


def cli_entry(argv: Sequence[str] | None = None) -> int:
    """Entry point for the console script."""
    try:
        main.main(args=argv, standalone_mode=False)
    except click.ClickException as exc:
        exc.show()
        return exc.exit_code
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(cli_entry())
