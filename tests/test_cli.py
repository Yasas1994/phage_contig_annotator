"""Unit tests for pca.cli."""

from pathlib import Path

from click.testing import CliRunner

from pca import cli, databases


def test_main_without_command_prints_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert result.exit_code == 0
    assert "usage:" in result.output.lower()


def test_run_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["run", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output
    assert "--translation-table" in result.output
    assert "--force" in result.output


def test_download_db_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["download-db", "--help"])
    assert result.exit_code == 0
    assert "--path" in result.output


def test_run_requires_input_and_output(tmp_path) -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["run"])
    assert result.exit_code != 0
    assert "Missing option" in result.output


def test_run_help_includes_gene_caller() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["run", "--help"])
    assert result.exit_code == 0
    assert "--gene-caller" in result.output
    assert "phanotate-rs" in result.output


def _make_minimal_db(db_dir: Path) -> None:
    db_dir.mkdir(parents=True, exist_ok=True)
    phrogs_dir = db_dir / "phrogs"
    (phrogs_dir / "hmms").mkdir(parents=True)
    (phrogs_dir / "hmms" / "test.hmm").write_text("dummy")
    (phrogs_dir / "metadata").mkdir(parents=True)
    (phrogs_dir / "metadata" / "PHROG_annot_v4.csv").write_text(
        "#phrog,Annotation,Category,color\n"
    )


def test_run_rejects_phanotate_when_not_installed(tmp_path: Path, monkeypatch) -> None:
    db_dir = tmp_path / "db"
    _make_minimal_db(db_dir)
    input_fna = tmp_path / "input.fna"
    input_fna.write_text(">contig1\n" + "ATG" * 20 + "TAA\n")

    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "run",
            "-i", str(input_fna),
            "-o", str(tmp_path / "out"),
            "--db", str(db_dir),
            "--gene-caller", "phanotate",
            "--skip-trna",
            "--skip-trf",
        ],
    )
    assert result.exit_code != 0
    assert "PHANOTATE is not on PATH" in result.output


def test_run_rejects_phanotate_rs_when_not_installed(tmp_path: Path, monkeypatch) -> None:
    db_dir = tmp_path / "db"
    _make_minimal_db(db_dir)
    input_fna = tmp_path / "input.fna"
    input_fna.write_text(">contig1\n" + "ATG" * 20 + "TAA\n")

    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "run",
            "-i", str(input_fna),
            "-o", str(tmp_path / "out"),
            "--db", str(db_dir),
            "--gene-caller", "phanotate-rs",
            "--skip-trna",
            "--skip-trf",
        ],
    )
    assert result.exit_code != 0
    assert "phanotate-rs is not on PATH" in result.output


def test_run_writes_phanotate_rs_gene_caller_to_config(tmp_path: Path, monkeypatch) -> None:
    import os
    import subprocess

    db_dir = tmp_path / "db"
    _make_minimal_db(db_dir)
    input_fna = tmp_path / "input.fna"
    input_fna.write_text(">contig1\n" + "ATG" * 20 + "TAA\n")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "phanotate-rs").write_text("#!/usr/bin/env python\nprint(1)\n")
    (bin_dir / "phanotate-rs").chmod(0o755)
    monkeypatch.setenv(
        "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
    )

    captured: list[dict] = []

    def fake_run(cmd, **kwargs):
        captured.append({"cmd": cmd, "kwargs": kwargs})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("pca.cli.subprocess.run", fake_run)

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "run",
            "-i", str(input_fna),
            "-o", str(tmp_path / "out"),
            "--db", str(db_dir),
            "--gene-caller", "phanotate-rs",
            "--skip-trna",
            "--skip-trf",
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured

    import yaml
    config_path = tmp_path / "out" / "config.yaml"
    assert config_path.is_file()
    config = yaml.safe_load(config_path.read_text())
    assert config["gene_caller"] == "phanotate-rs"


def test_run_writes_gene_caller_to_config(tmp_path: Path, monkeypatch) -> None:
    import os
    import subprocess

    db_dir = tmp_path / "db"
    _make_minimal_db(db_dir)
    input_fna = tmp_path / "input.fna"
    input_fna.write_text(">contig1\n" + "ATG" * 20 + "TAA\n")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "phanotate.py").write_text("#!/usr/bin/env python\nprint(1)\n")
    (bin_dir / "phanotate.py").chmod(0o755)
    monkeypatch.setenv(
        "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
    )

    captured: list[dict] = []

    def fake_run(cmd, **kwargs):
        captured.append({"cmd": cmd, "kwargs": kwargs})
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("pca.cli.subprocess.run", fake_run)

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "run",
            "-i", str(input_fna),
            "-o", str(tmp_path / "out"),
            "--db", str(db_dir),
            "--gene-caller", "phanotate",
            "--skip-trna",
            "--skip-trf",
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured

    import yaml
    config_path = tmp_path / "out" / "config.yaml"
    assert config_path.is_file()
    config = yaml.safe_load(config_path.read_text())
    assert config["gene_caller"] == "phanotate"


def test_run_converts_genbank_input_to_fasta(tmp_path: Path, monkeypatch) -> None:
    import os
    import subprocess

    from Bio import SeqIO
    from Bio.Seq import Seq

    db_dir = tmp_path / "db"
    _make_minimal_db(db_dir)

    input_gb = tmp_path / "input.fa"
    record = SeqIO.SeqRecord(
        Seq("ATG" * 20 + "TAA"),
        id="contig1",
        description="test phage",
        annotations={"molecule_type": "DNA"},
    )
    SeqIO.write(record, input_gb, "genbank")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "phanotate.py").write_text("#!/usr/bin/env python\nprint(1)\n")
    (bin_dir / "phanotate.py").chmod(0o755)
    monkeypatch.setenv(
        "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
    )

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("pca.cli.subprocess.run", fake_run)

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "run",
            "-i", str(input_gb),
            "-o", str(tmp_path / "out"),
            "--db", str(db_dir),
            "--gene-caller", "phanotate",
            "--skip-trna",
            "--skip-trf",
        ],
    )
    assert result.exit_code == 0, result.output

    import yaml
    config_path = tmp_path / "out" / "config.yaml"
    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    assert config["input"].endswith("input.fasta")
    fasta_path = Path(config["input"])
    assert fasta_path.is_file()
    assert ">contig1" in fasta_path.read_text()


def test_discover_databases_finds_flat_layout(tmp_path) -> None:
    (tmp_path / "test.hmm").write_text("dummy")
    (tmp_path / "PHROG_annot_v4.csv").write_text("dummy")
    discovered = databases.discover_databases(str(tmp_path))
    assert "phrogs" in discovered
    assert discovered["phrogs"].has_hmms
    assert discovered["phrogs"].has_metadata


def test_discover_databases_finds_hierarchical_layout(tmp_path) -> None:
    (tmp_path / "phrogs" / "hmms").mkdir(parents=True)
    (tmp_path / "phrogs" / "hmms" / "model.hmm").write_text("dummy")
    (tmp_path / "phrogs" / "metadata").mkdir(parents=True)
    (tmp_path / "phrogs" / "metadata" / "PHROG_annot_v4.csv").write_text("dummy")
    discovered = databases.discover_databases(str(tmp_path))
    assert "phrogs" in discovered
    assert discovered["phrogs"].has_hmms
    assert discovered["phrogs"].has_metadata


def test_discover_databases_finds_packaged_layout(tmp_path) -> None:
    (tmp_path / "hmmerdb_test").mkdir()
    (tmp_path / "hmmerdb_test" / "model.hmm").write_text("dummy")
    (tmp_path / "meta").mkdir()
    (tmp_path / "meta" / "PHROG_annot_v4.csv").write_text("dummy")
    discovered = databases.discover_databases(str(tmp_path))
    assert "hmmerdb" in discovered
    assert discovered["hmmerdb"].has_hmms
    assert discovered["hmmerdb"].has_metadata


def test_normalize_database_dir_organizes_flat_files(tmp_path) -> None:
    (tmp_path / "model.hmm").write_text("dummy")
    (tmp_path / "PHROG_annot_v4.csv").write_text("dummy")
    databases.normalize_database_dir(str(tmp_path))
    assert (tmp_path / "phrogs" / "hmms" / "model.hmm").is_file()
    assert (tmp_path / "phrogs" / "metadata" / "PHROG_annot_v4.csv").is_file()


def test_normalize_database_dir_organizes_multiple_databases(tmp_path) -> None:
    (tmp_path / "phrogs").mkdir()
    (tmp_path / "phrogs" / "9.hmm").write_text("dummy")
    (tmp_path / "phrogs" / "PHROG_annot_v4.csv").write_text("dummy")
    (tmp_path / "apis").mkdir()
    (tmp_path / "apis" / "dbAPIS.hmm").write_text("dummy")
    (tmp_path / "apis" / "dbAPIS.hmm.h3f").write_text("dummy")

    databases.normalize_database_dir(str(tmp_path))

    assert (tmp_path / "phrogs" / "hmms" / "9.hmm").is_file()
    assert (tmp_path / "phrogs" / "metadata" / "PHROG_annot_v4.csv").is_file()
    assert (tmp_path / "apis" / "hmms" / "dbAPIS.hmm").is_file()
    assert (tmp_path / "apis" / "hmms" / "dbAPIS.hmm.h3f").is_file()

    discovered = databases.discover_databases(str(tmp_path))
    assert "phrogs" in discovered
    assert "apis" in discovered
    assert discovered["phrogs"].has_hmms
    assert discovered["phrogs"].has_metadata
    assert discovered["apis"].has_hmms


    (tmp_path / "model.hmm").write_text("dummy")
    (tmp_path / "PHROG_annot_v4.csv").write_text("dummy")
    databases.normalize_database_dir(str(tmp_path))
    assert (tmp_path / "phrogs" / "hmms" / "model.hmm").is_file()
    assert (tmp_path / "phrogs" / "metadata" / "PHROG_annot_v4.csv").is_file()


def test_utils_reorganize_db_dry_run_lists_databases(tmp_path: Path) -> None:
    (tmp_path / "phrogs").mkdir()
    (tmp_path / "phrogs" / "9.hmm").write_text("dummy")
    (tmp_path / "phrogs" / "PHROG_annot_v4.csv").write_text("dummy")
    (tmp_path / "apis").mkdir()
    (tmp_path / "apis" / "dbAPIS.hmm").write_text("dummy")

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["utils", "reorganize-db", "--path", str(tmp_path), "--dry-run"],
    )
    assert result.exit_code == 0, result.output
    assert "phrogs" in result.output
    assert "apis" in result.output


def test_utils_reorganize_db_reorganizes_multiple_databases(tmp_path: Path) -> None:
    (tmp_path / "phrogs").mkdir()
    (tmp_path / "phrogs" / "9.hmm").write_text("dummy")
    (tmp_path / "phrogs" / "PHROG_annot_v4.csv").write_text("dummy")
    (tmp_path / "apis").mkdir()
    (tmp_path / "apis" / "dbAPIS.hmm").write_text("dummy")

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["utils", "reorganize-db", "--path", str(tmp_path)],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "phrogs" / "hmms" / "9.hmm").is_file()
    assert (tmp_path / "phrogs" / "metadata" / "PHROG_annot_v4.csv").is_file()
    assert (tmp_path / "apis" / "hmms" / "dbAPIS.hmm").is_file()


    runner = CliRunner()
    result = runner.invoke(cli.main, ["utils", "--help"])
    assert result.exit_code == 0
    assert "report" in result.output


def test_utils_report_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["utils", "report", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output
    assert "--output" in result.output
    assert "--fasta" in result.output


def test_utils_report_converts_genbank_to_html(tmp_path: Path) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord

    runner = CliRunner()

    gbk = tmp_path / "genome.gbk"
    record = SeqRecord(
        Seq("ATG" * 200 + "TAA"),
        id="contig1",
        annotations={"molecule_type": "DNA"},
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 300, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["gene_1"], "product": ["major head protein"]},
        )
    )
    SeqIO.write(record, gbk, "genbank")

    out_html = tmp_path / "report.html"
    result = runner.invoke(cli.main, ["utils", "report", "-i", str(gbk), "-o", str(out_html)])
    assert result.exit_code == 0, result.output
    assert out_html.is_file()
    content = out_html.read_text()
    assert "contig1" in content
    assert "major head protein" in content


def test_utils_report_rejects_single_html_for_multiple_records(tmp_path: Path) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    runner = CliRunner()

    gbk = tmp_path / "genome.gbk"
    records = [
        SeqRecord(Seq("ATG" * 100), id="c1", annotations={"molecule_type": "DNA"}),
        SeqRecord(Seq("ATG" * 100), id="c2", annotations={"molecule_type": "DNA"}),
    ]
    SeqIO.write(records, gbk, "genbank")

    out_html = tmp_path / "report.html"
    result = runner.invoke(cli.main, ["utils", "report", "-i", str(gbk), "-o", str(out_html)])
    assert result.exit_code != 0
    assert "specify an output directory" in result.output
