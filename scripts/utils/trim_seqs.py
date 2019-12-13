import argparse

from utils.bio import read_bio_seqs, write_bio_seqs


def trim_seqs(seqs, fr):
    return {s_id: s[int(fr*len(s)):int((1-fr)*len(s))]
            for s_id, s in seqs.items()}
