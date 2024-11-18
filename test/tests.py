import random
import subprocess
import time
import quicksect
from quicksect import Interval
from ncls import NCLS
import cgranges as cr
import superintervals

from superintervals import IntervalSet as superIntervalSet

import numpy as np
import pandas as pd
import sys


def make_random_bedtools(srt, l, n):
    random.seed(0)
    with open("chr1.genome", "w") as f:
        # f.write(f"chr1\t250000000")
        f.write(f"chr1\t10000000")  # 10 Mb

    subprocess.run(f"bedtools random -g chr1.genome -l 100 -n 1000000 -seed 1 > a.bed", shell=True)
    subprocess.run(f"bedtools random -g chr1.genome -l 50 -n 1000000 -seed 2 >> a.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l 1000000 -n 100 -seed 3 >> a.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l {l} -n {n} -seed 1 > a.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l {l*2} -n {n} -seed 2 >> a.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l {l*4} -n {n} -seed 3 >> a.bed", shell=True)

    subprocess.run(f"bedtools random -g chr1.genome -l {l} -n {n} -seed 4 > b.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l {l*2} -n {n} -seed 5 >> b.bed", shell=True)
    # subprocess.run(f"bedtools random -g chr1.genome -l {l*4} -n {n} -seed 6 >> b.bed", shell=True)

    intervals = []
    with open("a.bed", "r") as b:
        for line in b:
            l = line.split("\t")
            intervals.append((int(l[1]), int(l[2])))
    queries = []
    with open("b.bed", "r") as b:
        for line in b:
            l = line.split("\t")
            queries.append((int(l[1]), int(l[2])))
    if srt:
        queries.sort()
    intervals.sort()
    # subprocess.run("rm a.bed b.bed chr1.genome", shell=True)
    print("Made test intervals")
    return np.array(intervals), np.array(queries)


def make_worst_case(tower_size):
    # intervals, queries = make_worst_case(t)
    # res = run_tools(intervals, queries, srt)
    intervals = []
    s = 200_000
    e = 200_001
    for i in range(tower_size):
        intervals.append((s, e))
        s -= 100
        e += 100
    intervals.sort()
    queries = []
    for i in range(1):
        queries.append((i + 1_000_000, i + 1_000_000 + 1000))
    return np.array(intervals), np.array(queries)


def load_intervals(intervals_path, queries_path):
    queries = []
    intervals = []
    with open(intervals_path, "r") as f:
        for line in f:
            l = line.split("\t")
            intervals.append( (int(l[1]), int(l[2])) )
    with open(queries_path, "r") as f:
        for line in f:
            l = line.split("\t")
            queries.append( (int(l[1]), int(l[2])) )
    return np.array(intervals), np.array(queries)


def to_micro(t0):
    return int((time.time() - t0)*1000000)


def run_tools(intervals, queries):

    # superintervals
    t0 = time.time()
    sitv = superIntervalSet()
    for s, e in intervals:
        sitv.add(s, e, 0)
    sitv.index()
    print(f"SuperIntervals-py,{to_micro(t0)},",end='')
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = sitv.find_overlaps(start, end)
        v += len(a)
    print(f'{to_micro(t0)},{v},',end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        v += sitv.count_overlaps(start, end)
    print(f'{to_micro(t0)},{v}')

    # quicksect
    tree = quicksect.IntervalTree()
    for s, e in intervals:
        tree.add(s, e)
    print(f"Quicksect,{to_micro(t0)},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = tree.find(Interval(start, end))
        v += len(a)
    print(f'{to_micro(t0)},{v}')

    # cgranges
    cg = cr.cgranges()
    for s, e in intervals:
        cg.add("1", s, e+1, 0)
    cg.index()
    print(f"Cgranges-py,{to_micro(t0)},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(cg.overlap("1", start, end + 1))
        v += len(a)
    print(f'{to_micro(t0)},{v}')

    # ncls
    starts = pd.Series(intervals[:, 0])
    ends = pd.Series(intervals[:, 1])
    treencls = NCLS(starts, ends, starts)
    print(f"NCLS-py,{to_micro(t0)},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(treencls.find_overlap(start - 1, end + 1))
        v += len(a)
    print(f'{to_micro(t0)},{v}')


if __name__ == "__main__":

    intervals = sys.argv[1]
    queries = sys.argv[2]
    run_tools(*load_intervals(intervals, queries))
