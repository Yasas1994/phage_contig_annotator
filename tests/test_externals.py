"""Unit tests for pca.externals wrappers."""

from __future__ import annotations

import os
from pathlib import Path

from pca import externals, io


class TestSplitFasta:
    def test_split_fasta_creates_single_chunk_for_few_records(self, tmp_path: Path) -> None:
        fasta = tmp_path / "input.fna"
        fasta.write_text(
            ">c1\nACGT\n> c2\nTGCA\n> c3\nAAAA\n"
        )
        chunks = io.split_fasta(fasta, tmp_path / "chunks", 10)
        # n_chunks is capped at the number of records
        assert len(chunks) == 3
        all_text = "".join(p.read_text() for p in chunks)
        assert ">c1" in all_text
        assert "> c2" in all_text
        assert "> c3" in all_text

    def test_split_fasta_respects_requested_chunks(self, tmp_path: Path) -> None:
        fasta = tmp_path / "input.fna"
        fasta.write_text(
            ">c1\nACGT\n> c2\nTGCA\n> c3\nAAAA\n> c4\nCCCC\n> c5\nGGGG\n"
        )
        chunks = io.split_fasta(fasta, tmp_path / "chunks", 2)
        assert len(chunks) == 2
        assert {p.name for p in chunks} == {"chunk_0.fasta", "chunk_1.fasta"}
        # Sequential splitting: chunk 0 has c1, c2, c3; chunk 1 has c4, c5
        assert ">c1" in chunks[0].read_text()
        assert "> c5" in chunks[1].read_text()

    def test_split_fasta_empty_input(self, tmp_path: Path) -> None:
        fasta = tmp_path / "empty.fna"
        fasta.write_text("")
        chunks = io.split_fasta(fasta, tmp_path / "chunks", 4)
        assert len(chunks) == 1
        assert chunks[0].read_text() == ""


class TestTRF:
    @staticmethod
    def _make_fake_trf(bin_dir: Path) -> None:
        """Write a fake ``trf`` executable that mimics TRF output naming."""
        script = bin_dir / "trf"
        nl = chr(10)
        lines = [
            "#!/usr/bin/env python",
            "import sys, os",
            "in_file = sys.argv[1]",
            "params = sys.argv[2:9]",  # 7 numeric params
            "flags = sys.argv[9:]",  # -d -h etc.
            "out_file = in_file + '.' + '.'.join(params) + '.dat'",
            # Parse simple FASTA to emit one block per sequence header
            "contigs = []",
            "with open(in_file) as fh:",
            "    for line in fh:",
            "        line = line.strip()",
            "        if line.startswith('>'):",
            "            contigs.append(line[1:].split()[0])",
            "with open(out_file, 'w') as fh:",
            "    for c in contigs:",
            "        fh.write('Sequence: ' + c + '\\n')",
            "        fh.write('Parameters: 2 7 7 80 10 50 2000\\n')",
            "        fh.write('1 10 3 2.0 5 80 0 20 15 15 15 15 1.5 ATG ATGATG\\n')",
        ]
        script.write_text(nl.join(lines))
        script.chmod(0o755)

    def test_run_trf_single_chunk(self, tmp_path: Path, monkeypatch) -> None:
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        self._make_fake_trf(bin_dir)
        monkeypatch.setenv(
            "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        )

        fasta = tmp_path / "contigs.fna"
        fasta.write_text(
            ">c1\nACGTACGTACGTACGTACGTACGT\n"
            ">c2\nTGCATGCATGCATGCATGCATGCA\n"
        )

        out_prefix = str(tmp_path / "trf" / "trf")
        assert externals.run_trf(out_prefix, str(fasta), threads=1)

        out_dat = Path(f"{out_prefix}.dat")
        assert out_dat.is_file()
        text = out_dat.read_text()
        assert "Sequence: c1" in text
        assert "Sequence: c2" in text

    def test_run_trf_parallel_chunks(self, tmp_path: Path, monkeypatch) -> None:
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        self._make_fake_trf(bin_dir)
        monkeypatch.setenv(
            "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        )

        fasta = tmp_path / "contigs.fna"
        fasta.write_text(
            ">c1\nACGTACGTACGTACGTACGTACGT\n"
            ">c2\nTGCATGCATGCATGCATGCATGCA\n"
            ">c3\nAAAATTTTCCCCGGGGAAAA\n"
            ">c4\nGGGGAAAATTTTCCCC\n"
        )

        out_prefix = str(tmp_path / "trf" / "trf")
        assert externals.run_trf(out_prefix, str(fasta), threads=4)

        out_dat = Path(f"{out_prefix}.dat")
        assert out_dat.is_file()
        text = out_dat.read_text()
        assert text.count("Sequence:") == 4
        assert "Sequence: c1" in text
        assert "Sequence: c2" in text
        assert "Sequence: c3" in text
        assert "Sequence: c4" in text
        # Temporary chunk directory should have been removed
        assert not (tmp_path / "trf" / "trf_chunks").exists()

    def test_run_trf_creates_empty_dat_for_no_repeats(self, tmp_path: Path, monkeypatch) -> None:
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()

        # Fake trf that succeeds but writes no output file
        script = bin_dir / "trf"
        script.write_text("#!/usr/bin/env python\nimport sys\nprint('ok')\n")
        script.chmod(0o755)
        monkeypatch.setenv(
            "PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        )

        fasta = tmp_path / "contigs.fna"
        fasta.write_text(">c1\nACGTACGTACGTACGTACGTACGT\n")

        out_prefix = str(tmp_path / "trf" / "trf")
        assert externals.run_trf(out_prefix, str(fasta), threads=1)
        out_dat = Path(f"{out_prefix}.dat")
        assert out_dat.is_file()
        assert out_dat.stat().st_size == 0
