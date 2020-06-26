# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from sequence_graph.db_graph import DeBruijnGraph
from sequence_graph.idb_graph import get_db
from utils.various import fst_iterable

logger = logging.getLogger("centroFlye.sequence_graph.path_graph")


class PathDeBruijnGraph(DeBruijnGraph):
    def __init__(self, raw_pathdb, paths):
        self.__dict__.update(raw_pathdb.__dict__)
        self.paths = paths

    @classmethod
    def from_db(cls, db, string_set, k, neutral_symbs=None,
                assembly=None, outdir=None):
        if neutral_symbs is None:
            neutral_symbs = set()
        logger.info('Constructing a path graph')
        epaths = db.map_strings(string_set, only_unique_paths=True,
                                neutral_symbs=neutral_symbs)
        paths = {s_id: db.get_path(path, e_st=e_st, e_en=e_en)
                 for s_id, (path, e_st, e_en) in epaths.items()}
        logger.info(f'Extracted {len(paths)} paths')

        # TODO: rename. mode assembly is for min mult of a kmer == 1
        raw_pathdb, _ = get_db(paths, k, outdir,
                               assembly=assembly,
                               mode='assembly')
        pathdb = cls(raw_pathdb, paths)
        return pathdb

    @classmethod
    def from_mono_db(cls, db, monostring_set, k, assembly=None, outdir=None):
        monostring = fst_iterable(monostring_set.values())
        neutral_symbs = set([monostring.gap_symb])
        return cls.from_db(db=db,
                           string_set=monostring_set,
                           k=k,
                           neutral_symbs=neutral_symbs,
                           assembly=assembly,
                           outdir=outdir)
