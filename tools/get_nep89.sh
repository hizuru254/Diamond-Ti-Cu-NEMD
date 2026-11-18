#!/usr/bin/env bash
# Fetch the NEP89 model (14 MB) into potentials/.  SHA256 is recorded for the lock.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p potentials
if [ ! -f potentials/nep89.txt ]; then
  curl -sL -o potentials/nep89.txt \
    "https://raw.githubusercontent.com/brucefan1983/GPUMD/master/potentials/nep/nep89_20250409/nep89_20250409.txt"
fi
sha256sum potentials/nep89.txt
