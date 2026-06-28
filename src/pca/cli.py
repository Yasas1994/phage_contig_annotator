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

from pca import validation
from pca.logutils import get_logger


def _package_dir() -> Path:
    """Return the directory containing this package."""
    return Path(__file__).parent.resolve()


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


def _discover_databases(db_dir: str) -> tuple[dict[str, str], str | None]:
    """Return discovered HMM databases and the PHROG metadata path."""
    hmmerdb: dict[str, str] = {}
    meta_path: str | None = None

    for path in Path(db_dir).glob("*"):
        path_str = str(path)
        if "hmmerdb" in path_str:
            db_name = path.name.split("_")[0]
            hmmerdb[db_name] = path_str
        if "meta" in path_str:
            candidates = list(path.glob("phrog_annot_v4.tsv"))
            if candidates:
                meta_path = str(candidates[0])

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

    with tarfile.open(tar_file_path, "r:gz") as tar:
        for member in tar.getmembers():
            member.name = os.path.basename(member.name) or member.name
        tar.extractall(path=path)

    tar_file_path.unlink()
    click.echo(f"Database downloaded and extracted to {path}")
    checkpoint.touch()
    return True


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
    "--dry-run",
    is_flag=True,
    help="Show Snakemake execution plan without running.",
)
@click.pass_context
def runall(
    ctx: click.Context,
    input_path: Path,
    output_dir: Path,
    input_type: str,
    skip_trna: bool,
    db_dir: Path | None,
    cpus: int,
    plot_formats: tuple[str, ...],
    dry_run: bool,
) -> None:
    """Run the full annotation pipeline."""
    logger = get_logger(quiet=ctx.obj["quiet"])

    config = _load_config()
    db_dir = db_dir or Path(config["databases"]["dbroot"])
    if not db_dir.is_dir():
        raise click.UsageError(f"Database directory does not exist: {db_dir}")

    hmmerdb, meta_path = _discover_databases(str(db_dir))
    if not hmmerdb:
        raise click.UsageError(f"No HMM databases found in {db_dir}")
    if meta_path is None:
        raise click.UsageError(f"PHROG metadata not found in {db_dir}")

    if not skip_trna:
        validation.check_executables(["tRNAscan-SE"])

    output_dir.mkdir(parents=True, exist_ok=True)

    workflow_config = {
        "input": str(input_path.resolve()),
        "input_type": input_type,
        "output_dir": str(output_dir.resolve()),
        "db_dir": str(Path(next(iter(hmmerdb.values()))).resolve()),
        "meta_path": str(Path(meta_path).resolve()),
        "run_trna": not skip_trna,
        "plot_formats": list(plot_formats),
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

    env = os.environ.copy()
    src_dir = str(_package_dir().parent)
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{src_dir}{os.pathsep}{existing_pythonpath}" if existing_pythonpath else src_dir

    cmd.extend(["--envvars", "PYTHONPATH"])

    logger.info("running snakemake: %s", " ".join(cmd))
    result = subprocess.run(cmd, env=env)
    sys.exit(result.returncode)


@main.command()
@click.option(
    "-p",
    "--path",
    "path",
    default=None,
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    help="Path to store the database.",
)
@click.pass_context
def download_db(ctx: click.Context, path: Path | None) -> None:
    """Download the annotation database and update local config."""
    path = path or _package_dir() / "databases"
    path.mkdir(parents=True, exist_ok=True)

    if not download_dbs(path=str(path)):
        raise click.ClickException("Database download failed.")

    config = _load_config()
    config.set("databases", "dbroot", str(path))
    config.set("databases", "hmmdb", str(path / "hmmdb"))
    config.set("databases", "meta", str(path / "meta"))
    _write_config(config)
    click.echo(f"Updated config.ini with database path: {path}")


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
