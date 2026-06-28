#!/usr/bin/env bash
# One-line installer for phage_contig_annotator on macOS and Linux.
#
# Usage:
#   bash install.sh
#
# This script requires Conda (or Miniforge/Miniconda/Anaconda) to be installed.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$REPO_DIR/environment.yml"
ENV_NAME="phage_contig_annot"

# OS check
OS=$(uname -s)
case "$OS" in
    Linux|Darwin)
        ;;
    *)
        echo "ERROR: Unsupported operating system '$OS'. This installer supports Linux and macOS only." >&2
        exit 1
        ;;
esac

# Conda check
if ! command -v conda &>/dev/null; then
    echo "ERROR: Conda not found. Please install Miniforge, Miniconda, or Anaconda first:" >&2
    echo "  https://conda.io/projects/conda/en/latest/user-guide/install/index.html" >&2
    exit 1
fi

echo "==> Creating/updating conda environment '$ENV_NAME' from $ENV_FILE ..."
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    conda env update -f "$ENV_FILE" -n "$ENV_NAME"
else
    conda env create -f "$ENV_FILE"
fi

echo "==> Installing phage_contig_annotator in editable mode ..."
cd "$REPO_DIR"
conda run -n "$ENV_NAME" pip install -e ".[dev]"

echo "==> Downloading annotation database ..."
conda run -n "$ENV_NAME" --no-capture-output phage_contig_annotator download-db

echo ""
echo "==> Installation complete."
echo "Activate the environment with:"
echo "  conda activate $ENV_NAME"
echo "Then run:"
echo "  phage_contig_annotator run --input <fasta> --output <dir> --cpus <n>"
