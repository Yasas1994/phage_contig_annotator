"""Unit tests for pca.cli."""

from pathlib import Path

from click.testing import CliRunner

from pca import cli


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


def test_discover_databases_finds_flat_layout(tmp_path) -> None:
    (tmp_path / "test.hmm").write_text("dummy")
    (tmp_path / "PHROG_annot_v4.csv").write_text("dummy")
    hmmerdb, meta_path = cli._discover_databases(str(tmp_path))
    assert meta_path is not None
    assert Path(meta_path).name == "PHROG_annot_v4.csv"
    assert "hmmerdb" in hmmerdb


def test_discover_databases_finds_packaged_layout(tmp_path) -> None:
    (tmp_path / "hmmerdb_test").mkdir()
    (tmp_path / "hmmerdb_test" / "model.hmm").write_text("dummy")
    (tmp_path / "meta").mkdir()
    (tmp_path / "meta" / "PHROG_annot_v4.csv").write_text("dummy")
    hmmerdb, meta_path = cli._discover_databases(str(tmp_path))
    assert meta_path is not None
    assert "hmmerdb" in hmmerdb


def test_normalize_database_dir_organizes_flat_files(tmp_path) -> None:
    (tmp_path / "model.hmm").write_text("dummy")
    (tmp_path / "PHROG_annot_v4.csv").write_text("dummy")
    cli._normalize_database_dir(str(tmp_path))
    assert (tmp_path / "hmmerdb" / "model.hmm").is_file()
    assert (tmp_path / "meta" / "PHROG_annot_v4.csv").is_file()
