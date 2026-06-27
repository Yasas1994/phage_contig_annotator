"""Unit tests for pca.cli."""

from click.testing import CliRunner

from pca import cli


def test_main_without_command_prints_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, [])
    assert result.exit_code == 0
    assert "usage:" in result.output.lower()


def test_runall_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["runall", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output


def test_download_db_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["download-db", "--help"])
    assert result.exit_code == 0
    assert "--path" in result.output


def test_runall_requires_input_and_output(tmp_path) -> None:
    runner = CliRunner()
    result = runner.invoke(cli.main, ["runall"])
    assert result.exit_code != 0
    assert "Missing option" in result.output
