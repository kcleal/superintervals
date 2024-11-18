#!/bin/bash
#
# Overview
# ---------
# SuperIntervals (SI) was compared with:
# 1. Coitrees (Rust: https://github.com/dcjones/coitrees)
# 2. Implicit Interval Tree (C++: https://github.com/lh3/cgranges)
# 3. Interval Tree (C++: https://github.com/ekg/intervaltree)
# 4. Nested Containment List (C: https://github.com/pyranges/ncls/tree/master/ncls/src)
#
# Datasets:
# 1. Random regions generated using bedtools
# 2. RNA-seq reads and annotations from cgranges repository
# 3. ONT reads from sample PAO33946 (chr1 + chrM)
# 4. Paired-end reads from sample DB53, NCBI BioProject PRJNA417592, (chr1 + chrM)
# 5. ucsc genes for hg19
#
# Note, programs only assess chr1 bed records - other chromosomes are ignored. For 'chrM' records,
# the M was replaced with 1 using sed.
#
# Data were assessed in position sorted and random order. Datasets can be found on the Releases page.
#
# For Coitrees, a benchmark program is available at:
#     https://github.com/kcleal/superintervals/tree/main/test/3rd-party/coitrees/examples
# This will need to be build before hand

# Running:
# --------
# Run this script from the top-level superintervals repo
#     git clone https://github.com/kcleal/superintervals
#     cd superintervals
#
# Fetch data using:
#     wget https://github.com/kcleal/superintervals/releases/download/v0.2.0/data.tar.gz
#     tar -xvf data.tar.gz
#
# Run this script once:
#     bash test/run_tools.sh
#
# To run a few times and store results:
#     mkdir -p results
#     for r in {1..5}; do bash run_tools.sh 2>&1 | tee results/run${r}.txt done


run_libs=test/run-cpp-libs
rust_coitrees=test/3rd-party/coitrees/target/release/examples/bed-intersect
rust_superintervals=target/release/examples/bed-intersect-si


run_benchmarks() {
    local ref="$1"
    local query="$2"
    local label1="$3"
    local label2="$4"

    echo "$label1, $label2"
    $run_libs "$ref" "$query"
    $rust_superintervals "$ref" "$query"
    $rust_coitrees "$ref" "$query"
    $rust_coitrees -s "$ref" "$query"
    echo
}

# Prepare and run benchmark sets
prepare_and_run() {
    local ref="$1"
    local query="$2"
    local name1="$3"
    local name2="$4"

    # SORTED
    grep -P '^chr1\t' "$ref" | sort -k1,1 -k2,2n > ref
    grep -P '^chr1\t' "$query" | sort -k1,1 -k2,2n > query
    echo "SORTED"
    run_benchmarks ref query "$name1" "$name2"

    # SHUFFLED
    grep -P '^chr1\t' "$ref" | shuf > ref
    grep -P '^chr1\t' "$query" | shuf > query
    echo "SHUFFLED"
    run_benchmarks ref query "$name1" "$name2"

    echo "----------------------------------"
}

# random regions from bedtools random, baseline
prepare_and_run data/l1000_n1000.b.bed data/l1000_n1000.b.bed "rand-b" "rand-a"

# chrM reads from DB53 (renamed as chr1)
prepare_and_run data/DB53.chrM_reads_as_chr1.bed data/DB53.chrM_reads_as_chr1.bed "mito-b" "mito-a"
prepare_and_run data/PAO33946.chrM_reads_as_chr1.bed data/PAO33946.chrM_reads_as_chr1.bed "mito-lr-b" "mito-lr-a"

# RNA anno from cgranges releases page
prepare_and_run data/ex-rna.bed data/ex-anno.bed "rna" "anno"

# genes vs reads
prepare_and_run data/ucsc.hg19.genes.bed data/chr1.DB53.bed "genes" "DB53 reads"

# long reads vs short reads
prepare_and_run data/PAO33946.chr1.bed data/chr1.DB53.bed "ONT reads" "DB53 reads"

