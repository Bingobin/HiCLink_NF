#!/bin/bash
set -euo pipefail

usage() {
    printf '%s\n' \
        'Usage: bash rebuild_cooler.sh SAMPLE_ID VALID_PAIRS CHROM_SIZES [OUTDIR]' \
        'Defaults: OUTDIR=results/04_matrices, genome=hg38, CPUs=16, base=5000 bp.' \
        'Overrides: HICLINK_GENOME, HICLINK_CPUS, HICLINK_BASE_RESOLUTION,' \
        '           HICLINK_RESOLUTIONS, HICLINK_IGNORE_DIAGS.'
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    usage
    exit 0
fi
if (( $# < 3 || $# > 4 )); then
    usage >&2
    exit 2
fi

sample_id=$1
pairs=$2
chrom_sizes=$3
outdir=${4:-results/04_matrices}
genome=${HICLINK_GENOME:-hg38}
cpus=${HICLINK_CPUS:-16}
base_resolution=${HICLINK_BASE_RESOLUTION:-5000}
resolutions=${HICLINK_RESOLUTIONS:-5000,10000,25000,50000,100000,250000,500000,1000000}
ignore_diags=${HICLINK_IGNORE_DIAGS:-2}

[[ $sample_id =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || die 'Invalid sample ID.'
[[ $cpus =~ ^[1-9][0-9]*$ ]] || die 'HICLINK_CPUS must be a positive integer.'
[[ $base_resolution =~ ^[1-9][0-9]*$ ]] || die 'HICLINK_BASE_RESOLUTION must be a positive integer.'
[[ $ignore_diags =~ ^(0|[1-9][0-9]*)$ ]] || die 'HICLINK_IGNORE_DIAGS must be a non-negative integer.'
[[ $resolutions =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]] || die 'HICLINK_RESOLUTIONS must be comma-separated positive integers.'
IFS=',' read -r -a zoom_levels <<< "$resolutions"
for resolution in "${zoom_levels[@]}"; do
    (( resolution % base_resolution == 0 )) || die "Resolution $resolution is not a multiple of $base_resolution."
done
[[ ",$resolutions," == *",$base_resolution,"* ]] || die 'HICLINK_RESOLUTIONS must include the base resolution.'

command -v cooler >/dev/null 2>&1 || die 'cooler is not in PATH; activate the hiclink-nf Conda environment first.'
[[ -r $pairs && -s $pairs ]] || die "Missing or empty pairs file: $pairs"
[[ -r $chrom_sizes && -s $chrom_sizes ]] || die "Missing or empty chromosome sizes file: $chrom_sizes"

mkdir -p "$outdir"
base_cool="$outdir/$sample_id.$base_resolution.cool"
multi_cool="$outdir/$sample_id.mcool"
for output in "$base_cool" "$multi_cool"; do
    [[ ! -e $output && ! -L $output ]] || die "Refusing to overwrite $output; check it or choose a separate OUTDIR."
done

# Publish only after all three commands finish. Keep partial files for diagnosis.
build_dir=$(mktemp -d "$outdir/.$sample_id.build.XXXXXX")
trap 'build_status=$?; if (( build_status != 0 )); then printf "Build failed; intermediate files retained in %s\n" "$build_dir" >&2; fi' EXIT
staged_cool="$build_dir/$sample_id.$base_resolution.cool"
staged_mcool="$build_dir/$sample_id.mcool"

cooler cload pairs \
    --assembly "$genome" \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    "$chrom_sizes:$base_resolution" \
    "$pairs" \
    "$staged_cool"

cooler balance \
    --nproc "$cpus" \
    --ignore-diags "$ignore_diags" \
    "$staged_cool"

cooler zoomify \
    --balance \
    --balance-args "--nproc $cpus --ignore-diags $ignore_diags" \
    --nproc "$cpus" \
    --resolutions "$resolutions" \
    -o "$staged_mcool" \
    "$staged_cool"

cooler ls "$staged_mcool"
mv -n "$staged_cool" "$base_cool"
mv -n "$staged_mcool" "$multi_cool"
[[ ! -e $staged_cool && ! -e $staged_mcool ]] || die 'An output appeared during the build; existing files were not overwritten.'
rmdir "$build_dir"
printf 'Created:\n  %s\n  %s\n' "$base_cool" "$multi_cool"
