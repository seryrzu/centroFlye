# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from bisect import bisect_left


def dict_map(f, d):
    return {k: f(v) for k, v in d.items()}


def dict_map_name(f, d):
    return {k: f(k, v) for k, v in d.items()}


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0, myList[0]
    if pos == len(myList):
        return len(myList)-1, myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return pos, after
    else:
        return pos-1, before


# Taken from https://stackoverflow.com/a/4665027
def find_all_nonoverlap(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def find_all_overlap(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += 1  # use start += 1 to find overlapping matches


# Taken from https://stackoverflow.com/a/52177077
def chunks2(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


def get_kmers(fn):
    kmers = []
    with open(fn) as f:
        for line in f:
            kmers.append(line.strip())
    return set(kmers)


def list2str(lst, sep=' '):
    return sep.join(str(e) for e in lst)


def listEls2str(lst):
    return [str(e) for e in lst]


# from https://stackoverflow.com/a/40927266
def weighted_random_by_dct(dct):
    rand_val = np.random.random()
    total = 0
    for k, v in dct.items():
        total += v
        if rand_val <= total:
            return k
    assert False, 'unreachable'

def fst_iterable(iterable):
    return next(iter(iterable))
