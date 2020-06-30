# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import edlib

logger = logging.getLogger("centroFlye.mapping.mapping")


def map_queries(queries, targets,
                max_nloc_target,
                max_ntarget_locs,
                neutral_symbs,
                max_dist,
                mode='HW'):

    def map_query(query, target):
        if neutral_symbs is None or len(neutral_symbs) == 0:
            add_matches = None
        else:
            add_matches = []
            alphabet = set(query)
            for neutral_symb in neutral_symbs:
                for c in alphabet:
                    add_matches.append((c, neutral_symb))
        align = edlib.align(query,
                            target,
                            mode=mode,
                            task='locations',
                            k=max_dist,
                            additionalEqualities=add_matches)
        locs = align['locations']
        for i, (s, e) in enumerate(locs):
            locs[i] = (s, e+1)
        return locs

    all_locations = {i: {} for i in range(len(targets))}
    for q_id, query in queries.items():
        q_locs = {}
        for i, target in enumerate(targets):
            locs = map_query(query, target)
            if len(locs) > 0:
                q_locs[i] = locs
        if len(q_locs) <= max_ntarget_locs:
            for i, locs in q_locs.items():
                assert len(locs) > 0
                if len(locs) <= max_nloc_target:
                    all_locations[i][q_id] = locs
    return all_locations
