#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: build_deterministic_source.sh SOURCE OUTPUT_DIR" >&2
  exit 2
fi

source_dir=$(realpath "$1")
output_dir=$(realpath "$2")
mkdir -p "$output_dir"
scratch=$(mktemp -d /tmp/scwbGenePrograms-build-XXXXXX)
trap 'rm -rf -- "$scratch"' EXIT

(
  cd "$scratch"
  TZ=UTC R CMD build --no-build-vignettes --no-manual --compression=none \
    "$source_dir" >/dev/null
  mkdir unpacked
  tar -xf scwbGenePrograms_0.1.0.tar -C unpacked
)
sed -i '/^Packaged:/d' "$scratch/unpacked/scwbGenePrograms/DESCRIPTION"
tar --sort=name --mtime='@1' --owner=0 --group=0 --numeric-owner \
  --format=gnu -C "$scratch/unpacked" -cf - scwbGenePrograms | \
  gzip -n -9 > "$output_dir/scwbGenePrograms_0.1.0.tar.gz"
sha256sum "$output_dir/scwbGenePrograms_0.1.0.tar.gz"
