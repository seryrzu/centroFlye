# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
from copy import deepcopy
import logging

import numpy as np

from monomers.monostring import MonoString
from utils.bio import compress_homopolymer

logger = logging.getLogger("centroFlye.monomers.monostring_correction")


def get_stats(monostrings, return_stats=False):
    def get_ngap_symbols(monostrings, compr_hmp=False):
        cnt = 0
        for monostring in monostrings.values():
            string = monostring.string
            if compr_hmp:
                string = compress_homopolymer(string, return_list=True)
            char_counter = Counter(string)
            cnt += char_counter[monostring.gap_symb]
        return cnt

    stats = {}
    monostrings_lens = [len(monostr) for monostr in monostrings.values()]
    stats['ntranslations'] = len(monostrings_lens)
    stats['min_len'] = np.min(monostrings_lens)
    stats['max_len'] = np.max(monostrings_lens)
    stats['mean_len'] = np.mean(monostrings_lens)
    stats['tot_len'] = np.sum(monostrings_lens)
    stats['ngaps'] = get_ngap_symbols(monostrings)
    stats['pgaps'] = stats['ngaps'] / stats['tot_len']
    stats['ngap_runs'] = get_ngap_symbols(monostrings, compr_hmp=True)

    logger.info(f'Number of translations: {stats["ntranslations"]}')
    logger.info(f'Min length = {stats["min_len"]}')
    logger.info(f'Mean length = {stats["mean_len"]}')
    logger.info(f'Max length = {stats["max_len"]}')
    logger.info(f'Total length = {stats["tot_len"]}')

    logger.info(f'#(%) Gap symbols = {stats["ngaps"]} ({stats["pgaps"]})')
    logger.info(f'#Gap runs = {stats["ngap_runs"]}')
    if return_stats:
        return stats


# def filter_lowercaserich_strings(monostrings, max_lowercase=0.1):
#     def filter_lowercaserich_string(monostring):
#     filtered_strings = {}


def correct_monostrings(raw_monostrings, inplace=True):
    if not inplace:
        raw_monostrings = deepcopy(raw_monostrings)
    get_stats(raw_monostrings)
    monostrings = raw_monostrings
    return monostrings
