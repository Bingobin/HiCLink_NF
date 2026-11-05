#!/usr/bin/env bash
set -euo pipefail

# Build the current official FitHiChIP scripts on top of the pinned official runtime.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-$(mktemp -d "${TMPDIR:-/tmp}/fithichip-docker.XXXXXX")}"
commit="0ea1ac21be870908c672316ffbb630189dc6fae2"
mkdir -p "$build_dir/FitHiChIP"

curl -fL --retry 3 "https://codeload.github.com/ay-lab/FitHiChIP/tar.gz/$commit" -o "$build_dir/source.tar.gz"
tar -xzf "$build_dir/source.tar.gz" --strip-components=1 -C "$build_dir/FitHiChIP"
curl -fL --retry 3 'https://repo.anaconda.com/miniconda/Miniconda3-py310_25.1.1-2-Linux-x86_64.sh' \
    -o "$build_dir/Miniconda3-py310_25.1.1-2-Linux-x86_64.sh"
curl -fL --retry 3 'https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.14.8.tar.gz' \
    -o "$build_dir/data.table_1.14.8.tar.gz"
curl -fL --retry 3 'https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedToBigBed' -o "$build_dir/bedToBigBed"
curl -fL --retry 3 'https://genome.ucsc.edu/goldenPath/help/examples/interact/interact.as' -o "$build_dir/interact.as"

printf '%s  %s\n' \
    7f298109ab95b5436632973a04189a125282cc948f1dd1b03fa9cb6c71443915 "$build_dir/Miniconda3-py310_25.1.1-2-Linux-x86_64.sh" \
    14b2ce5367df9c9bb58f373555066f5dcb629c156149b5565de36d69557139fd "$build_dir/data.table_1.14.8.tar.gz" \
    | shasum -a 256 --check
cp "$repo_dir/docker/fithichip/Dockerfile" "$build_dir/Dockerfile"
docker build --platform linux/amd64 -t hiclink/fithichip:11-0ea1ac2 "$build_dir"
printf 'FitHiChIP image installed. Build files retained at: %s\n' "$build_dir"
