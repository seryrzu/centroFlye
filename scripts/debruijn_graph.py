from collections import defaultdict

import networkx as nx
import numpy as np


class DeBruijnGraph:
    def __init__(self, k, max_uniq_cov=30):
        self.graph = nx.MultiDiGraph()
        self.k = k
        self.node_mapping = {}
        self.rev_node_mapping = {}
        self.max_uniq_cov = max_uniq_cov

    def add_kmer(self, kmer, coverage=1):
        prefix, suffix = kmer[:-1], kmer[1:]

        if prefix in self.node_mapping:
            prefix_node_ind = self.node_mapping[prefix]
        else:
            prefix_node_ind = len(self.node_mapping)
            self.node_mapping[prefix] = prefix_node_ind
            self.rev_node_mapping[prefix_node_ind] = prefix

        if suffix in self.node_mapping:
            suffix_node_ind = self.node_mapping[suffix]
        else:
            suffix_node_ind = len(self.node_mapping)
            self.node_mapping[suffix] = suffix_node_ind
            self.rev_node_mapping[suffix_node_ind] = suffix

        self.graph.add_edge(prefix_node_ind, suffix_node_ind,
                            edge_kmer=kmer,
                            length=1,
                            coverages=[coverage],
                            label=f'len=1\ncov={coverage}',
                            color='black' if coverage < self.max_uniq_cov else 'red')

    def add_kmers(self, kmers, coverage=None):
        for kmer in kmers:
            if coverage is None:
                self.add_kmer(kmer)
            else:
                self.add_kmer(kmer, coverage=coverage[kmer])

    def index_edges(self, min_k=2):
        all_index = {}
        for k in range(min_k, self.k+1):
            index = defaultdict(list)
            for e_ind, edge in enumerate(self.graph.edges(keys=True)):
                edge_seq = self.graph.get_edge_data(*edge)['edge_kmer']
                for i in range(len(edge_seq)-k+1):
                    kmer = edge_seq[i:i+k]
                    index[kmer].append((e_ind, i))
            index_unique = {kmer: pos[0]
                            for kmer, pos in index.items()
                            if len(pos) == 1}
            all_index[k] = index_unique
        return all_index

    def collapse_nonbranching_paths(self):
        def node_on_nonbranching_path(graph, node):
            return nx.number_of_nodes(graph) > 1 \
                and graph.in_degree(node) == 1 \
                and graph.out_degree(node) == 1

        for node in list(self.graph.nodes()):
            if node_on_nonbranching_path(self.graph, node):
                in_edge = list(self.graph.in_edges(node, keys=True))[0]
                out_edge = list(self.graph.out_edges(node, keys=True))[0]
                # in_edge_color = self.graph.edges[in_edge]['color']
                # out_edge_color = self.graph.edges[out_edge]['color']
                in_edge_kmer = self.graph.edges[in_edge]['edge_kmer']
                out_edge_kmer = self.graph.edges[out_edge]['edge_kmer']
                in_edge_cov = self.graph.edges[in_edge]['coverages']
                out_edge_cov = self.graph.edges[out_edge]['coverages']

                in_node = in_edge[0]
                out_node = out_edge[1]

                new_kmer = in_edge_kmer + \
                    out_edge_kmer[-(len(out_edge_kmer)-self.k+1):]
                new_coverages = in_edge_cov + out_edge_cov
                new_coverages.sort()
                new_edge_len = len(new_coverages)
                new_edge_med_cov = np.median(new_coverages)
                self.graph.add_edge(in_node, out_node,
                                    edge_kmer=new_kmer,
                                    coverages=new_coverages,
                                    length=new_edge_len,
                                    label=f'len={new_edge_len}\ncov={new_edge_med_cov}',
                                    color='black' if new_edge_med_cov < self.max_uniq_cov else 'red')
                self.graph.remove_node(node)

    def get_edges(self):
        self.collapse_nonbranching_paths()
        contigs, coverages = [], []
        for u, v, keys in self.graph.edges(data=True):
            contigs.append(keys['edge_kmer'])
            coverages.append(np.median(keys['coverages']))
        return contigs, coverages

    def get_path(self, list_edges):
        db_edges = list(self.graph.edges())
        path = ''
        path += self.graph.get_edge_data(*list_edges[0])['edge_kmer']
        for edge in list_edges[1:]:
            new_edge_seq = self.graph.get_edge_data(*edge)['edge_kmer']
            # print(len(path[-(self.k-1):]), path[-(self.k-1):])
            # print(len(new_edge_seq[:self.k-1]), new_edge_seq[:self.k-1])
            assert path[-(self.k-1):] == new_edge_seq[:self.k-1]
            path += new_edge_seq[self.k-1:]

        # process cycled path
        if list_edges[0][0] == list_edges[-1][1]:
            path = path[:-(self.k-1)]
            doubled_path = path*2
            min_shift, min_path = 0, path
            for i in range(1, len(path)):
                cur_path = doubled_path[i:i+len(path)]
                if cur_path < min_path:
                    min_shift = i
                    min_path = cur_path
            path = min_path
        return path

    def get_contigs(self):
        def get_valid_outpaths_from_edge(edge, taken_edges=set()):
            out_node = edge[1]
            out_degree_out_node = self.graph.out_degree(out_node)
            paths = [[edge]]
            if out_degree_out_node == 1:
                out_edge = list(self.graph.out_edges(out_node, keys=True))[0]
                if out_edge not in taken_edges:
                    passed_taken_edges = taken_edges.copy()
                    passed_taken_edges.add(edge)
                    paths = get_valid_outpaths_from_edge(out_edge,
                                                         taken_edges=passed_taken_edges)
                    paths = paths.copy()
                    for i in range(len(paths)):
                        paths[i].append(edge)
            return paths

        def get_valid_inpaths_from_edge(edge, taken_edges=set()):
            in_node = edge[0]
            in_degree_in_node = self.graph.in_degree(in_node)

            paths = [[edge]]
            if in_degree_in_node == 1:
                in_edge = list(self.graph.in_edges(in_node, keys=True))[0]
                if in_edge not in taken_edges:
                    passed_taken_edges = taken_edges.copy()
                    passed_taken_edges.add(edge)
                    paths = get_valid_inpaths_from_edge(in_edge,
                                                        taken_edges=passed_taken_edges)
                    paths = paths.copy()
                    for i in range(len(paths)):
                        paths[i].append(edge)
            return paths

        # Getting valid in- and out-paths for each edge
        edge_outpaths, edge_inpaths = {}, {}
        for edge in self.graph.edges(keys=True):
            edge_outpaths[edge] = get_valid_outpaths_from_edge(edge)
            edge_inpaths[edge] = get_valid_inpaths_from_edge(edge)
        for edge in self.graph.edges(keys=True):
            edge_inpaths[edge] = [tuple(path) for path in edge_inpaths[edge]]
            edge_outpaths[edge] = [tuple(path[::-1]) for path in edge_outpaths[edge]]

        # Unite all such paths and delete duplicates
        all_paths = []
        for edge in self.graph.edges(keys=True):
            all_paths += edge_inpaths[edge] + edge_outpaths[edge]
        all_paths = list(set(all_paths))

        # Select only paths that are not a subpath of other paths
        selected_paths = []
        # doubled_paths = [path + path for path in all_paths]
        for path1 in all_paths:
            has_duplicate = False
            for path2 in all_paths:
                if path1 == path2:
                    continue
                for i in range(len(path2)-len(path1)+1):
                    if path1 == path2[i:i+len(path1)]:
                        has_duplicate = True
                        break
                if has_duplicate:
                    break
            if not has_duplicate:
                selected_paths.append(path1)

        # check that all edges of that graph are incorporated in selected paths
        # edges_in_selected_paths = []
        # for path in selected_paths:
        #     edges_in_selected_paths += path
        # assert len(set(self.graph.edges(keys=True)) - set(edges_in_selected_paths)) == 0

        contigs = []
        for path in selected_paths:
            contigs.append(self.get_path(path))
        contigs = list(set(contigs))
        return contigs, selected_paths
    

def iterative_graph(monomer_strings, min_k, max_k, outdir, min_mult=5, step=1, starting_graph=None):
    smart_makedirs(outdir)
    
    dbs, all_contigs = {}, {}
    input_strings = monomer_strings.copy()
    
    if starting_graph is not None:
        contigs, contig_paths = starting_graph.get_contigs()
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{min_k}_i{i}_j{j}'] = contigs[i]

    for k in range(min_k, max_k+1, step):
        print(f'\nk={k}')
        # print(len(input_strings))
        frequent_kmers, frequent_kmers_read_pos = get_frequent_kmers(input_strings, k=k, min_mult=min_mult)
        print(f'#frequent kmers = {len(frequent_kmers)}')
        db = DeBruijnGraph(k=k)
        db.add_kmers(frequent_kmers, coverage=frequent_kmers)

        # lens = sorted((len(contig), coverage) for contig, coverage in zip(contigs, coverages))[::-1]
        # long_edges = [x for x in lens if x[1] >= 50]
        
        db.collapse_nonbranching_paths()
        if nx.number_weakly_connected_components(db.graph) > 1:
            print(f'#cc = {nx.number_weakly_connected_components(db.graph)}')
            for cc in nx.weakly_connected_components(db.graph):
                print(len(cc))
            # break
        dbs[k] = db
        
        dot_file = os.path.join(outdir, f'db_k{k}.dot')
        pdf_file = os.path.join(outdir, f'db_k{k}.pdf')
        nx.drawing.nx_pydot.write_dot(db.graph, dot_file)
        !dot -Tpdf {dot_file} -o {pdf_file}
        
        
        contigs, contig_paths = db.get_contigs()
        all_contigs[k] = contigs

        input_strings = monomer_strings.copy()
        # print(len(input_strings))
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{k}_i{i}_j{j}'] = contigs[i]
        # print(len(input_strings))
    
    return all_contigs, dbs
