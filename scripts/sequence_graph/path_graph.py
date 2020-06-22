# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from sequence_graph.db_graph import DeBruijnGraph
from sequence_graph.idb_graph import get_db

logger = logging.getLogger("centroFlye.sequence_graph.path_graph")


class PathDeBruijnGraph(DeBruijnGraph):
    def __init__(self, raw_pathdb):
        self.__dict__.update(raw_pathdb.__dict__)

    @classmethod
    def from_db(cls, db, strings, k, outdir=None):
        logger.info('Constructing a path graph')
        mappings = db.map_strings(strings)
        paths = {s_id: db.get_path(mapping.epath)
                 for s_id, mapping in mappings.items()
                 if mapping is not None and mapping.valid}
        logger.info(f'Extracted {len(paths)} paths')

        # TODO: rename. mode assembly is for min mult of a kmer == 1
        raw_pathdb, _ = get_db(strings, k, outdir, mode='assembly')
        pathdb = cls(raw_pathdb)
        return pathdb
