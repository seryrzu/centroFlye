# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import edlib

from monomers.monostring import MonoInstance, MonoString

logger = logging.getLogger("centroFlye.submonomers.submonostring")


class SubmonoInstance(MonoInstance):
    # submonomers can be empty if it is unidentified
    def __init__(self, submonomers, dist2submonomers,
                 monomer, strand, seq_id, nucl_segment,
                 st, en, seq_len,
                 reliability):
        super().__init__(monomer, strand, seq_id, nucl_segment,
                         st, en, seq_len, reliability)
        self.submonomers = submonomers
        self.dist2submonomers = dist2submonomers

    @classmethod
    def from_monoinstance(cls, monoinstance, submonomers, dist2submonomers):
        return cls(submonomers=submonomers,
                   dist2submonomers=dist2submonomers,
                   monomer=monoinstance.monomer,
                   strand=monoinstance.strand,
                   seq_id=monoinstance.seq_id,
                   nucl_segment=monoinstance.nucl_segment,
                   st=monoinstance.st,
                   en=monoinstance.en,
                   seq_len=monoinstance.seq_len,
                   reliability=monoinstance.reliability)

    def is_unequivocal(self):
        return len(self.submonomers) == 1

    def is_unidentified(self):
        return len(self.submonomers) == 0

    def is_identified(self):
        identified = not self.is_unidentified()
        if identified:
            assert self.is_reliable
        return identified

    def is_equivocal(self):
        # return True if len(self.submonomers) > 1
        return self.is_identified() and not self.is_unequivocal()

    def make_unequivocal(self, submonoindex):
        submonomers_dict = {submonomer.submonoindex: submonomer
                            for submonomer in self.submonomers}
        assert submonoindex in submonomers_dict.keys()
        self.submonomers = [submonomers_dict[submonoindex]]


class SubmonoString(MonoString):
    submonogap_symb = '¿'    # if monomer is identifed, but no submonomer
    submonoequiv_symb = '*'  # if several submonomers are identified

    def __init__(self, seq_id, monoinstances, raw_monostring, nucl_sequence,
                 monomer_db, is_reversed,
                 submonoinstances, raw_submonostring, submonomer_db):
        super().__init__(seq_id=seq_id,
                         monoinstances=monoinstances,
                         raw_monostring=raw_monostring,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed)
        assert len(submonoinstances) == len(raw_submonostring)
        self.submonoinstances = submonoinstances
        self.raw_submonostring = raw_submonostring
        self.submonomer_db = submonomer_db

    @classmethod
    def from_monostring(cls, monostring, submonomer_db,
                        max_k=5, ext=10, max_diff=3):
        raw_submonostring = []
        submonoinstances = []

        monogap_symb = MonoString.gap_symb
        submonogap_symb = SubmonoString.submonogap_symb
        submonoequiv_symb = SubmonoString.submonoequiv_symb

        for i, mi in enumerate(monostring.monoinstances):
            assert len(submonoinstances) == len(raw_submonostring) == i
            if not mi.is_reliable():
                raw_submonostring.append(monogap_symb)
                sminst = SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                           submonomers={},
                                                           dist2submonomers={})
                submonoinstances.append(sminst)
                continue

            st, en = mi.st - ext, mi.en + ext
            nucl_segment = monostring.get_nucl_segment(st, en)
            monoindex = mi.get_monoindex()
            submonomers = \
                submonomer_db.get_submonomers_by_mono_index(monoindex)

            dist_smis = []
            for sm in submonomers:
                alignment = edlib.align(sm.seq, nucl_segment,
                                        k=max_k, mode='HW')
                dist = alignment['editDistance']
                if dist == -1:
                    continue
                dist_smis.append((dist, sm.submono_index))

            dist_smis.sort()
            if len(dist_smis) == 0:
                raw_submonostring.append(submonogap_symb)
                sminst = SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                           submonomers={},
                                                           dist2submonomers={})
                submonoinstances.append(sminst)
                continue

            distances, smindexes = list(zip(*dist_smis))
            best_dist, best_smindex = distances[0], smindexes[0]
            if len(distances) == 1:
                raw_submonostring.append(best_smindex)
                sm = submonomer_db.submonomers[best_smindex]
                submonomers = {best_smindex: sm}
                d2s = {best_smindex: best_dist}
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
                continue

            sec_best_dist, sec_best_smindex = distances[1], smindexes[1]
            # len(distances) > 1
            if best_dist == sec_best_dist == 0:
                best_sm = submonomer_db.submonomers[best_smindex]
                sec_best_sm = submonomer_db.submonomers[sec_best_smindex]
                logger.error(f'distance 0 to at least 2 segments:\n'
                             f'seq id = {monostring.seq_id}\n'
                             f'index in the string = {i}\n'
                             f'monoindex = {monoindex}\n'
                             f'nucl segment = {nucl_segment}\n'
                             f'best submonomer = {best_smindex}, '
                             f'{best_sm.seq}'
                             f'sec best submonomer = {sec_best_smindex}',
                             f'{sec_best_sm.seq}')
            if best_dist == 0 and sec_best_dist > 0:
                raw_submonostring.append(best_smindex)
                sm = submonomer_db.submonomers[best_smindex]
                submonomers = {best_smindex: sm}
                # distance to submonomer
                d2s = {best_smindex: best_dist}  # best_dist == 0
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
            else:
                submonomers = {}
                d2s = {}  # distance to submonomer
                for distance, smindex in zip(distances, smindexes):
                    if distance - best_dist > max_diff:
                        # distances is sorted
                        break
                    d2s[smindex] = distance
                    sm = submonomer_db.submonomers[smindex]
                    submonomers[smindex] = sm
                raw_submonostring.append(submonoequiv_symb)
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
        submonostring = cls(seq_id=monostring.seq_id,
                            monoinstances=monostring.monoinstances,
                            raw_monostring=monostring.raw_monostring,
                            nucl_sequence=monostring.nucl_sequence,
                            monomer_db=monostring.monomer_db,
                            is_reversed=monostring.is_reversed,
                            submonoinstances=submonoinstances,
                            raw_submonostring=raw_submonostring,
                            submonomer_db=submonomer_db)
        return submonostring
