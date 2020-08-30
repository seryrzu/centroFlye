# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import os

from utils.os_utils import smart_makedirs
from utils.various import fst_iterable

from mapping.mapping import map_queries, get_coverage


def map_monoreads_to_monoassembly(monoreads_set, monoassembly_set,
                                  outdir=None, plot_close=True, title=None,
                                  max_dist=0):
    monoassembly = fst_iterable(monoassembly_set.values())
    raw_monoassembly = monoassembly.raw_monostring
    gap_symb = monoassembly.gap_symb
    locations = map_queries(queries=monoreads_set,
                            targets=[raw_monoassembly],
                            max_nloc_target=1,
                            max_ntarget_locs=1,
                            neutral_symbs=set(gap_symb),
                            max_dist=max_dist)
    locations = locations[0]

    if outdir is None:
        return locations

    smart_makedirs(outdir)
    locations_fn = os.path.join(outdir, 'locations.tsv')
    with open(locations_fn, 'w') as f:
        print(f'r_id\tmono_st\tmono_en\tr_nucl_st\tr_nucl_en\t'
              f'a_nucl_st\ta_nucl_en', file=f)
        for r_id, loc in locations.items():
            assert len(loc) == 1
            s, e = loc[0]
            a_nucl_s = monoassembly.monoinstances[s].st
            a_nucl_e = monoassembly.monoinstances[e-1].en
            monoread = monoreads_set[r_id]
            r_nucl_s = monoread.monoinstances[0].st
            r_nucl_e = monoread.monoinstances[-1].en
            print(f'{r_id}\t{s}\t{e}\t{r_nucl_s}\t{r_nucl_e}\t'
                  f'{a_nucl_s}\t{a_nucl_e}', file=f)

    if title is None:
        title = 'Coverage of monoassembly with monoreads'
    get_coverage(locations=locations,
                 target_len=len(raw_monoassembly),
                 outdir=outdir,
                 title=title,
                 plot_close=plot_close)
    return locations


def map_paths_to_monoassembly(paths, monoassembly_set, outdir=None,
                              plot_close=True, title=None,
                              max_dist=0):
    monoassembly = fst_iterable(monoassembly_set.values())
    raw_monoassembly = monoassembly.raw_monostring
    locations = map_queries(queries=paths,
                            targets=[raw_monoassembly],
                            max_nloc_target=1,
                            max_ntarget_locs=1,
                            neutral_symbs=None,
                            max_dist=max_dist)
    locations = locations[0]

    if outdir is None:
        return locations

    smart_makedirs(outdir)
    locations_fn = os.path.join(outdir, 'locations.tsv')
    with open(locations_fn, 'w') as f:
        print(f'path_id\tmono_st\tmono_en', file=f)
        for r_id, loc in locations.items():
            assert len(loc) == 1
            s, e = loc[0]
            print(f'{r_id}\t{s}\t{e}', file=f)

    if title is None:
        title = 'Coverage of monoassembly with paths'
    get_coverage(locations=locations,
                 target_len=len(raw_monoassembly),
                 outdir=outdir,
                 title=title,
                 plot_close=plot_close)
    return locations
