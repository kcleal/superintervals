import random
import subprocess
from sortedintersect import IntervalSet
import time
import quicksect
from quicksect import Interval
from ncls import NCLS
import cgranges as cr
from superintervals import IntervalSet as superIntervalSet
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

random.seed(0)

def make_random_bedtools(shuffle, n=1_000_000, l=1000):

    with open("chr1.genome", "w") as f:
        f.write(f"chr1\t250000000")

    subprocess.run(f"bedtools random -g chr1.genome -l {l} -n {n} -seed 1 > a.bed", shell=True)
    subprocess.run(f"bedtools random -g chr1.genome -l {l*2} -n {n} -seed 2 >> a.bed", shell=True)
    subprocess.run(f"bedtools random -g chr1.genome -l {20000000} -n {10} -seed 3 >> a.bed", shell=True)

    subprocess.run(f"bedtools random -g chr1.genome -l {l} -n {n} -seed 4 > b.bed", shell=True)
    subprocess.run(f"bedtools random -g chr1.genome -l {l*2} -n {n} -seed 5 >> b.bed", shell=True)
    subprocess.run(f"bedtools random -g chr1.genome -l {20000000} -n {10} -seed 6 >> b.bed", shell=True)

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
    if not shuffle:
        queries.sort()
    intervals.sort()
    # subprocess.run("rm a.bed b.bed chr1.genome", shell=True)
    print("Made test intervals")
    return np.array(intervals), np.array(queries)


def load_intervals(shuffle):
    dirname = os.path.dirname(__file__)
    queries = []
    intervals = []
    with open(dirname + "/chr1_ucsc_genes.bed", "r") as f:
        for line in f:
            l = line.split("\t")
            intervals.append( (int(l[1]), int(l[2])) )
    with open(dirname + "/chr1_reads.bed", "r") as f:
        for line in f:
            l = line.split("\t")
            queries.append( (int(l[1]), int(l[2])) )
    if shuffle:
        random.shuffle(queries)
    else:
        queries.sort()
    intervals.sort()
    return np.array(intervals), np.array(queries)


def run_tools(intervals, queries, shuffled):

    res = []

    # superintervals
    sitv = superIntervalSet()
    for s, e in intervals:
        sitv.add(s, e, 0)
    sitv.index()

    # sortedintersect
    itv = IntervalSet(False)
    itv.add_from_iter(intervals)

    # quicksect
    tree = quicksect.IntervalTree()
    for s, e in intervals:
        tree.add(s, e)

    # cgranges
    cg = cr.cgranges()
    for s, e in intervals:
        cg.add("1", s, e+1, 0)
    cg.index()

    # ncls
    starts = pd.Series(intervals[:, 0])
    ends = pd.Series(intervals[:, 1])
    treencls = NCLS(starts, ends, starts)

    print("Interval sets constructed")
    print("superintervals-count")
    t0 = time.time()
    v = 0
    for start, end in queries:
        v += sitv.count_overlaps(start, end)
    res.append({'library': 'superintervals-count', 'time (s)': time.time() - t0, 'intersections': v, 'queries': len(queries)})
    print("superintervals-find")
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = sitv.find_overlaps(start, end)
        v += len(a)
    res.append({'library': 'superintervals-find', 'time (s)': time.time() - t0, 'intersections': v, 'queries': len(queries)})


    # t0 = time.time()
    # v = 0
    # for start, end in queries:
    #     a = itv.search_interval(start, end)
    #     v += len(a)
    # res.append({'library': 'sortedintersect', 'time (s)': time.time() - t0, 'intersections': v})
    print("quicksect")
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = tree.find(Interval(start, end))
        # print(a)
        v += len(a)
    res.append({'library': 'quicksect', 'time (s)': time.time() - t0, 'intersections': v, 'queries': len(queries)})

    print("cgranges")
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(cg.overlap("1", start, end+1))
        v += len(a)
    res.append({'library': 'cgranges', 'time (s)': time.time() - t0, 'intersections': v, 'queries': len(queries)})

    print("ncls")
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(treencls.find_overlap(start - 1, end + 1))
        v += len(a)
    res.append({'library': 'ncls', 'time (s)': time.time() - t0, 'intersections': v, 'queries': len(queries)})


    return pd.DataFrame.from_records(res)



dfs = []
for shuffle in (True,): #, False):
    print("Shuffled data", shuffle)
    print("random1")
    intervals, queries = make_random_bedtools(shuffle, n=1_000_000, l=100)
    res = run_tools(intervals, queries, shuffle)
    res["test"] = ["random1"] * len(res)
    res["shuffle"] = [shuffle] * len(res)
    dfs.append(res)
    # print("random2")
    # intervals, queries = make_random_bedtools(shuffle, n=1_000_000, l=100)
    # df2 = run_tools(intervals, queries, shuffle)
    # df2["test"] = ["random2"] * len(df2)
    # df2["shuffle"] = [shuffle] * len(df2)
    # dfs.append(df2)

    # print("reads+genes")
    # intervals, queries = load_intervals(shuffle)
    # df2 = run_tools(intervals, queries, shuffle)
    # df2["test"] = ["reads+genes"] * len(df2)
    # df2["shuffle"] = [shuffle] * len(df2)
    # dfs.append(df2)
    print()


df = pd.concat(dfs)
print(df.to_markdown(index=False))

sns.set_palette("Set2")
sns.catplot(kind='bar', data=df, x="shuffle", y="time (s)", hue="library", col="test", sharey=False)
plt.show()
# plt.savefig('benchmark.png')

