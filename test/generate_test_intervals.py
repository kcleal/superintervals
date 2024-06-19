"""
This script required bedtools to be available on your PATH.
Run this script from within this test directory
"""

from subprocess import run
import numpy as np
from scipy import stats
import random
import os


g = "genome.txt"
out = "bench_data"
seed = 0

if not os.path.exists(out):
    os.mkdir(out)


def generate_lognormal(mean, variance, size):
    # Calculate mu and sigma
    mu = np.log(mean ** 2 / np.sqrt(variance + mean ** 2))
    sigma = np.sqrt(np.log(1 + variance / mean ** 2))

    # Generate log-normal distributed random variables
    return stats.lognorm(sigma, loc=mean, scale=np.exp(mu)).rvs(size=size)


def existing_overlaps(target_queries):
    """All queries overlap at least one reference interval"""
    rand_sizes_q = generate_lognormal(1000, 5000, target_queries).astype(int)
    rand_sizes_r = generate_lognormal(1000, 5000, target_queries).astype(int)
    rand_pos = [random.randint(1, 250_000_000) for _ in range(target_queries)]
    with open(f"{out}/REF_query_overlaps_all_ref.{target_queries}.bed", "w") as out_ref, \
         open(f"{out}/QUERY_query_overlaps_all_ref.{target_queries}.bed", "w") as out_query:
        for sq, sr, pos in zip(rand_sizes_q, rand_sizes_r, rand_pos):
            end = pos + sq
            out_query.write(f"chr1\t{pos}\t{end}\n")
            pos_ref = random.randint(max(1, pos - sr - 1), end - 1)
            out_ref.write(f"chr1\t{pos_ref}\t{pos_ref + sr}\n")

    run(f"bedtools sort -i {out}/REF_query_overlaps_all_ref.{target_queries}.bed > {out}/REF_query_overlaps_all_ref.{target_queries}.srt.bed", shell=True)
    run(f"bedtools sort -i {out}/QUERY_query_overlaps_all_ref.{target_queries}.bed > {out}/QUERY_query_overlaps_all_ref.{target_queries}.srt.bed", shell=True)

    print(f"Target queries {target_queries}, unsorted, then sorted:")
    run(f"./run-cpp-libs {out}/REF_query_overlaps_all_ref.{target_queries}.bed {out}/QUERY_query_overlaps_all_ref.{target_queries}.bed", shell=True)
    run(f"./run-cpp-libs {out}/REF_query_overlaps_all_ref.{target_queries}.srt.bed {out}/QUERY_query_overlaps_all_ref.{target_queries}.srt.bed", shell=True)


existing_overlaps(1_000_000)


def non_existing_overlaps(target_queries):
    """All queries overlap zero reference intervals"""
    rand_sizes_q = [1] * target_queries*2
    rand_sizes_r = [1] * target_queries*2
    rand_pos = [random.randint(1, 250_000_000) for _ in range(target_queries)]
    with open(f"{out}/REF_query_overlaps_0_ref.{target_queries}.tmp.bed", "w") as out_ref, \
            open(f"{out}/QUERY_query_overlaps_0_ref.{target_queries}.bed", "w") as out_query:
        for sq, sr, pos in zip(rand_sizes_q, rand_sizes_r, rand_pos):
            end = pos + sq
            out_query.write(f"chr1\t{pos}\t{end}\n")
            pos_ref = end + 1
            out_ref.write(f"chr1\t{pos_ref}\t{pos_ref + sr}\n")

    run(f"bedtools subtract -a {out}/REF_query_overlaps_0_ref.{target_queries}.tmp.bed -b {out}/QUERY_query_overlaps_0_ref.{target_queries}.bed > {out}/REF_query_overlaps_0_ref.{target_queries}.bed", shell=True)
    run(f"rm {out}/REF_query_overlaps_0_ref.{target_queries}.tmp.bed", shell=True)

    run(f"bedtools sort -i {out}/REF_query_overlaps_0_ref.{target_queries}.bed > {out}/REF_query_overlaps_0_ref.{target_queries}.srt.bed",
        shell=True)
    run(f"bedtools sort -i {out}/QUERY_query_overlaps_0_ref.{target_queries}.bed > {out}/QUERY_query_overlaps_0_ref.{target_queries}.srt.bed",
        shell=True)

    print(f"Target queries {target_queries}, unsorted, then sorted:")
    run(f"./run-cpp-libs {out}/REF_query_overlaps_0_ref.{target_queries}.bed {out}/QUERY_query_overlaps_0_ref.{target_queries}.bed",
        shell=True)
    run(f"./run-cpp-libs {out}/REF_query_overlaps_0_ref.{target_queries}.srt.bed {out}/QUERY_query_overlaps_0_ref.{target_queries}.srt.bed",
        shell=True)

# non_existing_overlaps(1_000_000)