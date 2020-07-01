# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import os

import ahocorasick
import networkx as nx

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from mapping.mapping import map_queries
from utils.os_utils import smart_makedirs
from utils.various import fst_iterable


class AutomatonGraph:
    def __init__(self, gr, jumps, st_node):
        self.gr = gr
        self.jumps = jumps
        self.st_node = st_node

    @classmethod
    def from_automaton(cls, automaton):
        gr = nx.DiGraph()
        nodes, edges, jumps = automaton.dump()
        st_node = nodes[0][0]
        for node, marker in nodes:
            gr.add_node(node, marker=marker)
        gr.node[st_node]['len'] = 0
        for i, (s, c, e) in enumerate(edges):
            s_len = gr.nodes[s]['len']
            c = c.decode()
            gr.nodes[e]['len'] = s_len + 1
            gr.add_edge(s, e, char=c)
        jumps = {s: e for s, e in jumps}
        return cls(gr=gr, jumps=jumps, st_node=st_node)

    def jump_node(self, s, c):
        for s1, e in self.gr.out_edges(s):
            assert s == s1
            edge_c = self.gr.get_edge_data(s, e)['char']
            if c == edge_c:
                return e

    def get_substr_survival(self, raw_monoassembly):
        substr_surv = [0] * len(raw_monoassembly)
        node = self.st_node
        for i, c in enumerate(raw_monoassembly):
            e = None
            e = self.jump_node(s=node, c=c)
            if e is None:
                # can't extend, need to use the AC-links in trie
                length = self.gr.nodes[node]['len']
                substr_surv[i-length] = length
                while e is None and node != self.st_node:
                    node = self.jumps[node]
                    e = self.jump_node(s=node, c=c)
            if e is not None:
                node = e

        for i in range(1, len(raw_monoassembly)):
            substr_surv[i] = max(substr_surv[i],
                                 substr_surv[i-1]-1,
                                 0)
        return substr_surv


def _monolist2monostr(monolist):
    return ''.join(chr(i) if isinstance(i, int) else i for i in monolist)


def find_longest_surviving_substring(mapped_strings, monoassembly_set, outdir,
                                     plot_close=True):
    monoassembly = fst_iterable(monoassembly_set.values())
    raw_monoassembly = monoassembly.raw_monostring
    automaton = ahocorasick.Automaton()
    for s_id, string in mapped_strings.items():
        automaton.add_word(_monolist2monostr(string), (s_id, string))
    automaton.make_automaton()
    automaton_gr = AutomatonGraph.from_automaton(automaton)
    substr_surv = \
        automaton_gr.get_substr_survival(_monolist2monostr(raw_monoassembly))

    smart_makedirs(outdir)
    substr_surv_fn = os.path.join(outdir, 'substr_surv.tsv')
    with open(substr_surv_fn, 'w') as f:
        for i, length in enumerate(substr_surv):
            print(f'{i}\t{length}', file=f)

    substr_surv_pdf_fn = os.path.join(outdir, 'substr_surv.pdf')
    plt.plot(substr_surv)
    plt.xlabel('monoassembly')
    plt.ylabel('substring length')
    plt.title('Assembly substring survival')
    plt.grid()
    plt.savefig(substr_surv_pdf_fn, format='pdf')
    if plot_close:
        plt.close()


def find_longest_surviving_substring_paths(paths, monoassembly_set, outdir,
                                           plot_close=True):
    monoassembly = fst_iterable(monoassembly_set.values())
    raw_monoassembly = monoassembly.raw_monostring
    locations = map_queries(paths, [raw_monoassembly],
                            max_nloc_target=1, max_ntarget_locs=1,
                            neutral_symbs=None, max_dist=0)[0]
    mapped_paths = {r_id: paths[r_id] for r_id in locations}
    find_longest_surviving_substring(mapped_strings=mapped_paths,
                                     monoassembly_set=monoassembly_set,
                                     outdir=outdir,
                                     plot_close=True)
