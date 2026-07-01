"""Command-line interface for the phage contig annotator pipeline."""

from __future__ import annotations

import configparser
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Any, Sequence

import click
import requests
import tqdm

from pca import annotations, databases, genome_stats, validation
from pca.io import convert_to_fasta
from pca.logutils import get_logger


def _package_dir() -> Path:
    """Return the directory containing this package."""
    return Path(__file__).parent.resolve()


def _normalize_input_fasta(input_path: Path, output_dir: Path) -> Path:
    """Ensure the pipeline input is a FASTA file.

    GenBank (.gb, .gbk) and EMBL inputs are converted to a multi-FASTA file
    inside ``output_dir``. Already-FASTA inputs are returned unchanged.
    """
    from pca.io import detect_sequence_format

    fmt = detect_sequence_format(input_path)
    if fmt == "fasta":
        return input_path

    normalized = output_dir / "input.fasta"
    convert_to_fasta(input_path, normalized)
    return normalized


_DEFAULT_DB_DIR: Path = _package_dir().parents[1] / "databases"

# Optional extra databases searched with HMMER (hmmsearch for .hmm profiles,
# phmmer for protein FASTA databases). DefenseFinder is run as an external
# tool but exposed through the same interface.
_EXTRA_DB_NAMES = [
    "card",
    "vfdb",
    "anticrisprdb",
    "netflax",
    "dgr",
    "defensefinder",
]


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
        "3NbY2CPn74M5C9R/download/database.tar.gz"
    )
    os.makedirs(path, exist_ok=True)

    response = requests.get(url, stream=True, timeout=(30, 600))
    if response.status_code != 200:
        click.echo(
            f"Failed to download the database from {url}. "
            f"Status code: {response.status_code}",
            err=True,
        )
        return False

    tar_file_path = Path(path) / "database.tar.gz"
    total_size = int(response.headers.get("content-length", 0)) or None
    with tqdm.tqdm(
        total=total_size,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        desc="database.tar.gz",
        miniters=1,
        dynamic_ncols=True,
        bar_format="{desc}: {percentage:3.0f}%|{bar:5}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
    ) as pbar:
        with open(tar_file_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    file.write(chunk)
                    pbar.update(len(chunk))

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

    databases.normalize_database_dir(path)

    tar_file_path.unlink()
    click.echo(f"Database downloaded and extracted to {path}")
    checkpoint.touch()
    return True


def _primary_database(
    db_dir: Path, primary_name: str = "phrogs"
) -> databases.DatabaseBundle:
    """Discover databases and return the primary bundle.

    Raises a click.UsageError when the primary database is missing or lacks
    HMMs/metadata.
    """
    discovered = databases.discover_databases(str(db_dir))
    if not discovered:
        raise click.UsageError(
            f"No databases found in {db_dir}. "
            "Run 'phage_contig_annotator download-db' first."
        )

    primary = discovered.get(primary_name)
    if primary is None:
        available = ", ".join(sorted(discovered))
        raise click.UsageError(
            f"Primary database '{primary_name}' not found in {db_dir}. "
            f"Available databases: {available}"
        )

    if not primary.has_hmms:
        raise click.UsageError(
            f"Primary database '{primary_name}' has no HMM profiles in "
            f"{primary.hmms_path}"
        )
    if not primary.has_metadata:
        raise click.UsageError(
            f"Primary database '{primary_name}' has no metadata file."
        )

    return primary


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
    help="Path to input FASTA, GenBank (.gb, .gbk) or EMBL file with contigs or proteins.",
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
    "--gene-caller",
    "gene_caller",
    default="pyrodigal-gv",
    show_default=True,
    type=click.Choice(["pyrodigal-gv", "phanotate", "phanotate-rs"], case_sensitive=False),
    help="Gene prediction method for contigs (ignored for protein inputs).",
)
@click.option(
    "--skip-trna",
    is_flag=True,
    help="Skip tRNAscan-SE gene prediction.",
)
@click.option(
    "--skip-trf",
    is_flag=True,
    help="Skip Tandem Repeats Finder repeat detection.",
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
    "--extra-db-root",
    "extra_db_root",
    default=None,
    type=click.Path(file_okay=False, readable=True, path_type=Path),
    help="Root directory for optional extra databases. Defaults to <db>/extra.",
)
@click.option(
    "--extra-db",
    "extra_dbs",
    multiple=True,
    type=click.Choice(_EXTRA_DB_NAMES, case_sensitive=False),
    help="Enable an optional extra database. Specify multiple times.",
)
@click.option(
    "--extra-eval",
    "extra_evalue",
    default=1e-5,
    show_default=True,
    type=float,
    help="E-value threshold for optional extra database searches.",
)
@click.option(
    "--run-defensefinder",
    is_flag=True,
    help="Run DefenseFinder to detect anti-phage defense systems.",
)
@click.option(
    "--run-crisprcasfinder",
    is_flag=True,
    help="Run CRISPRCasFinder to detect CRISPR arrays and Cas genes.",
)
@click.option(
    "--run-minced",
    is_flag=True,
    help="Run MinCED to detect CRISPR arrays.",
)
@click.option(
    "--skip-multicopy",
    is_flag=True,
    help="Skip multi-copy detection on contigs (it runs by default).",
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
    gene_caller: str,
    skip_trna: bool,
    skip_trf: bool,
    db_dir: Path | None,
    extra_db_root: Path | None,
    extra_dbs: tuple[str, ...],
    extra_evalue: float,
    run_defensefinder: bool,
    run_crisprcasfinder: bool,
    run_minced: bool,
    skip_multicopy: bool,
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
    primary_name = config["databases"].get("primary", "phrogs")
    if not primary_name:
        primary_name = "phrogs"
    if not db_dir.is_dir():
        raise click.UsageError(f"Database directory does not exist: {db_dir}")

    primary = _primary_database(db_dir, primary_name=primary_name)

    if not skip_trna:
        validation.check_executables(["tRNAscan-SE"])
    run_trf = (not skip_trf) and input_type == "contigs"
    if run_trf:
        validation.check_executables(["trf"])

    if run_crisprcasfinder:
        if not shutil.which("CRISPRCasFinder.pl") and not shutil.which("CRISPRCasFinder"):
            raise click.UsageError(
                "CRISPRCasFinder was requested but CRISPRCasFinder.pl is not on PATH."
            )
    if run_minced:
        if not shutil.which("minced"):
            raise click.UsageError("MinCED was requested but minced is not on PATH.")

    if not skip_multicopy and input_type != "contigs":
        raise click.UsageError("Multi-copy detection is only supported for contig inputs.")

    gene_caller = gene_caller.lower()
    if input_type == "contigs":
        if gene_caller == "phanotate" and not (
            shutil.which("phanotate.py") or shutil.which("phanotate")
        ):
            raise click.UsageError(
                "Gene caller 'phanotate' was selected but PHANOTATE is not on PATH. "
                "Install it (`pip install phanotate`) or choose 'pyrodigal-gv'."
            )
        if gene_caller == "phanotate-rs" and not shutil.which("phanotate-rs"):
            raise click.UsageError(
                "Gene caller 'phanotate-rs' was selected but phanotate-rs is not on PATH. "
                "Install it (https://github.com/Yasas1994/PHANOTATE-rs) or choose 'pyrodigal-gv'."
            )

    extra_dbs_list = list(extra_dbs)
    if run_defensefinder and "defensefinder" not in extra_dbs_list:
        extra_dbs_list.append("defensefinder")
    if "defensefinder" in extra_dbs_list:
        run_defensefinder = True

    extra_db_root_path = extra_db_root or db_dir / "extra"
    if extra_dbs_list:
        for db in extra_dbs_list:
            if db == "defensefinder":
                validation.check_executables(["defense-finder"])
                continue
            db_path = extra_db_root_path / db
            if not db_path.exists():
                raise click.UsageError(
                    f"Extra database directory does not exist: {db_path}"
                )
            has_models = (
                any(db_path.glob("*.hmm"))
                or any(
                    p.suffix.lower() in {".faa", ".fasta", ".fa"}
                    for p in db_path.iterdir()
                )
            )
            if not has_models:
                raise click.UsageError(
                    f"Extra database {db_path} must contain .hmm files or a protein FASTA"
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = _normalize_input_fasta(input_path, output_dir)

    workflow_config: dict[str, Any] = {
        "input": str(input_path.resolve()),
        "input_type": input_type,
        "gene_caller": gene_caller,
        "output_dir": str(output_dir.resolve()),
        "db_dir": str(db_dir.resolve()),
        "primary_db": primary.name,
        "meta_path": str(primary.metadata_path.resolve()),
        "run_trna": not skip_trna,
        "run_trf": run_trf,
        "plot_formats": list(plot_formats),
        "translation_table": translation_table,
        "theme": theme,
        "extra_db_root": str(extra_db_root_path.resolve()) if extra_db_root else None,
        "extra_dbs": extra_dbs_list,
        "extra_evalue": extra_evalue,
        "run_defensefinder": run_defensefinder,
        "run_crisprcasfinder": run_crisprcasfinder,
        "run_minced": run_minced,
        "run_multicopy": not skip_multicopy,
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
    config.set("databases", "primary", "phrogs")
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
@click.option(
    "--stats",
    "stats_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Optional per-contig statistics TSV to display in the HTML report.",
)
@click.pass_context
def utils_report(
    ctx: click.Context,
    input_path: Path,
    output_path: Path,
    fasta_path: Path | None,
    meta_path: Path | None,
    theme: str,
    stats_tsv: Path | None,
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
            stats_tsv=stats_tsv,
        )
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    for out_path in out_paths:
        click.echo(f"Wrote {out_path}")
    logger.info("generated %d HTML report(s)", len(out_paths))


@utils.command("stats")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Path to input nucleotide FASTA file.",
)
@click.option(
    "-g",
    "--gff",
    "gff_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Path to GFF file with gene/CDS features.",
)
@click.option(
    "-t",
    "--trf",
    "trf_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Optional pre-computed TRF .dat file.",
)
@click.option(
    "--run-trf",
    is_flag=True,
    help="Run TRF with flags '-f -d -m' if no --trf file is provided.",
)
@click.option(
    "--skip-multicopy",
    is_flag=True,
    help="Skip multi-copy detection on the input FASTA (it runs by default).",
)
@click.option(
    "--multicopy-tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Optional pre-computed multi-copy detection TSV to merge into the stats table.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    default=None,
    type=click.Path(file_okay=True, dir_okay=False, writable=True, path_type=Path),
    help="Output TSV path (default: stdout).",
)
@click.pass_context
def utils_stats(
    ctx: click.Context,
    input_path: Path,
    gff_path: Path,
    trf_path: Path | None,
    run_trf: bool,
    skip_multicopy: bool,
    multicopy_tsv: Path | None,
    output_path: Path | None,
) -> None:
    """Compute per-contig genomic statistics (strand switching, repeats, copy number)."""
    logger = get_logger(quiet=ctx.obj["quiet"])

    trf_dat = trf_path
    if trf_dat is None and run_trf:
        from pca.externals import run_trf

        tmp_dir = Path(output_path).parent if output_path else Path.cwd()
        tmp_dir.mkdir(parents=True, exist_ok=True)
        trf_prefix = tmp_dir / f"{input_path.stem}_trf"
        if not run_trf(str(trf_prefix), str(input_path), flags=["-f", "-d", "-m"]):
            raise click.ClickException("TRF failed; cannot compute tandem-repeat stats.")
        trf_dat = Path(f"{trf_prefix}.dat")

    multicopy_path = multicopy_tsv
    if multicopy_path is None and not skip_multicopy:
        from pca.multicopy import detect_multicopy

        tmp_dir = Path(output_path).parent if output_path else Path.cwd()
        tmp_dir.mkdir(parents=True, exist_ok=True)
        multicopy_path = tmp_dir / f"{input_path.stem}_multicopy.tsv"
        df_mc = detect_multicopy(str(input_path))
        df_mc.to_csv(multicopy_path, sep="\t", index=False)
        logger.info("wrote multi-copy detection TSV to %s", multicopy_path)

    try:
        df = genome_stats.compute_genome_stats(
            input_path,
            gff_path,
            trf_path=trf_dat,
            multicopy_path=multicopy_path,
        )
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    if output_path is None:
        click.echo(df.to_csv(sep="\t", index=False))
    else:
        df.to_csv(output_path, sep="\t", index=False)
        click.echo(f"Wrote {output_path}")
    logger.info("computed stats for %d contig(s)", len(df))


@utils.command("reorganize-db")
@click.option(
    "-p",
    "--path",
    "db_dir",
    default=None,
    type=click.Path(exists=True, file_okay=False, writable=True, readable=True, path_type=Path),
    help="Path to the database directory to reorganize (default: bundled databases dir).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show the planned layout without moving files.",
)
@click.pass_context
def utils_reorganize_db(
    ctx: click.Context,
    db_dir: Path | None,
    dry_run: bool,
) -> None:
    """Reorganize an existing database directory into the hierarchical layout.

    Each top-level subdirectory becomes a database entry, with HMM profiles,
    metadata, DIAMOND and MMseqs files moved into standard ``hmms/``,
    ``metadata/``, ``diamond/`` and ``mmseqs/`` subdirectories.
    """
    logger = get_logger(quiet=ctx.obj["quiet"])
    db_dir = db_dir or _DEFAULT_DB_DIR

    before = databases.discover_databases(str(db_dir))
    if not before:
        raise click.ClickException(f"No databases found in {db_dir}.")

    if dry_run:
        click.echo(f"Planned layout for {db_dir}:")
        for name, bundle in sorted(before.items()):
            click.echo(f"  {name}/")
            if bundle.hmms_path is not None:
                click.echo(f"    hmms/ -> {bundle.hmms_path}")
            if bundle.metadata_path is not None:
                click.echo(f"    metadata/ -> {bundle.metadata_path}")
            if bundle.diamond_path is not None:
                click.echo(f"    diamond/ -> {bundle.diamond_path}")
            if bundle.mmseqs_path is not None:
                click.echo(f"    mmseqs/ -> {bundle.mmseqs_path}")
        return

    reorganized = databases.normalize_database_dir(str(db_dir))
    click.echo(f"Reorganized database directory: {db_dir}")
    for name, bundle in sorted(reorganized.items()):
        click.echo(f"  {name}: hmms={bundle.hmms_path}, metadata={bundle.metadata_path}")
    logger.info("reorganized database directory %s", db_dir)



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
