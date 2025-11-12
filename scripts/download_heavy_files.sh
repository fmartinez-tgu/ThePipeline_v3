#!/usr/bin/env bash
# Download helper for large files that were removed from the repo to save space.
#
# Usage:
#  - Option A (automated): set one or more of the following environment variables to direct-download URLs
#       MINOS_URL         -> Singularity image for Minos (minos_v0.12.5.img)
#       GATK_LOCAL_URL    -> gatk-package-4.2.5.0-local.jar URL
#       GATK_SPARK_URL    -> gatk-package-4.2.5.0-spark.jar URL
#       SNIPEFF_BIN_URL   -> snpEffectPredictor.bin URL (or use snpEff's internal download)
#    then run: ./scripts/download_heavy_files.sh
#
#  - Option B (manual): run the script with no env vars set. It will print exact paths and instructions
#    (including commands you can run manually after you have obtained the files via browser/portal).
#
# IMPORTANT:
#  - GATK requires you to accept Broad's license and download from the Broad website. Automated download
#    may not work for GATK. If you cannot download GATK automatically, please download it manually and
#    place the jar(s) under Programs/gatk-4.2.5.0/.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MINOS_DEST="$REPO_ROOT/Programs/Minos/minos_v0.12.5.img"
GATK_DIR="$REPO_ROOT/Programs/gatk-4.2.5.0"
GATK_LOCAL_DEST="$GATK_DIR/gatk-package-4.2.5.0-local.jar"
GATK_SPARK_DEST="$GATK_DIR/gatk-package-4.2.5.0-spark.jar"
SNIPEFF_BIN_DEST="$REPO_ROOT/Programs/snpEff/data/GRCh37.75/snpEffectPredictor.bin"

mkdir -p "$(dirname "$MINOS_DEST")" "${GATK_DIR}" "$(dirname "$SNIPEFF_BIN_DEST")"

download_if_var() {
  local url="$1"; shift
  local dest="$1"; shift
  if [ -n "${url:-}" ]; then
    echo "Downloading to $dest from $url ..."
    if command -v curl >/dev/null 2>&1; then
      curl -L --fail --progress-bar -o "$dest" "$url"
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$dest" "$url"
    else
      echo "Neither curl nor wget found; cannot download automatically. Please download $url manually and place it at $dest." >&2
      return 2
    fi
    echo "Download finished: $dest"
    return 0
  else
    return 1
  fi
}

echo
echo "=== ThePipeline3 heavy files download helper ==="
echo "Repository root: $REPO_ROOT"
echo

EXIT_CODE=0

echo "-- Minos Singularity image --"
if download_if_var "${MINOS_URL:-}" "$MINOS_DEST"; then
  echo "Minos image downloaded to $MINOS_DEST"
else
  echo "No MINOS_URL provided. To obtain the Singularity image, either:
  1) Download the Minos Singularity image from the Minos release page (GitHub Releases) if available,
     or
  2) Build the Singularity image locally from a Docker image with:
       singularity build $MINOS_DEST docker://<minos-docker-image>

Place the final file at: $MINOS_DEST
" | sed 's/^/    /'
fi

echo
echo "-- GATK jars --"
if download_if_var "${GATK_LOCAL_URL:-}" "$GATK_LOCAL_DEST"; then
  echo "Downloaded GATK local jar to $GATK_LOCAL_DEST"
else
  echo "No GATK_LOCAL_URL provided. GATK distributions require accepting Broad Institute's license.
  Please download the GATK jars manually from the Broad Institute (https://github.com/broadinstitute/gatk/releases or https://software.broadinstitute.org/gatk/download)
  After downloading, place the file(s) in: $GATK_DIR

  Expected file names:
    - gatk-package-4.2.5.0-local.jar -> $GATK_LOCAL_DEST
    - gatk-package-4.2.5.0-spark.jar -> $GATK_SPARK_DEST
" | sed 's/^/    /'
fi

if download_if_var "${GATK_SPARK_URL:-}" "$GATK_SPARK_DEST"; then
  echo "Downloaded GATK spark jar to $GATK_SPARK_DEST"
fi

echo
echo "-- snpEff data binary (GRCh37.75) --"
if download_if_var "${SNIPEFF_BIN_URL:-}" "$SNIPEFF_BIN_DEST"; then
  echo "Downloaded snpEff data binary to $SNIPEFF_BIN_DEST"
else
  echo "No SNIPEFF_BIN_URL provided. Two options to obtain snpEffectPredictor.bin for GRCh37.75:
  1) If you have snpEff installed (snpEff jar), run:
       java -jar snpEff.jar download -v GRCh37.75
     This will place the data under the snpEff data directory (commonly <snpEff_dir>/data/GRCh37.75/)

  2) Or download the dataset directly from the snpEff/SourceForge or dataset mirror and place
     the resulting 'snpEffectPredictor.bin' at:
       $SNIPEFF_BIN_DEST

After placing the file, ensure proper ownership and permissions.
" | sed 's/^/    /'
fi

echo
echo "Files present in Programs/ (summary):"
ls -lh "$REPO_ROOT/Programs/Minos" || true
ls -lh "$GATK_DIR" || true
ls -lh "$(dirname "$SNIPEFF_BIN_DEST")" || true

echo
echo "Done. If automatic downloads failed because you did not supply URLs, please download the files manually and re-run this script or set the environment variables:
  MINOS_URL, GATK_LOCAL_URL, GATK_SPARK_URL, SNIPEFF_BIN_URL

Examples (automated):
  MINOS_URL=https://example.com/minos_v0.12.5.img \
  GATK_LOCAL_URL=https://example.com/gatk-package-4.2.5.0-local.jar \
  GATK_SPARK_URL=https://example.com/gatk-package-4.2.5.0-spark.jar \
  SNIPEFF_BIN_URL=https://example.com/snpEffectPredictor.bin \
  ./scripts/download_heavy_files.sh

Notes:
  - GATK distributed jars often require manual agreement to licensing terms; automated downloads may not be permitted.
  - Minos Singularity images may be large; ensure sufficient disk space.

Exit code: $EXIT_CODE

exit $EXIT_CODE
