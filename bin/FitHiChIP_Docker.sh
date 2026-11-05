#!/usr/bin/env bash
set -euo pipefail

# Preserve the official -C interface; run all scientific dependencies inside Docker.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$script_dir/fithichip.py" run-call "$@"
