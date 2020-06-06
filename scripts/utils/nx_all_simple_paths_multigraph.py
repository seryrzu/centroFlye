# This is from
# https://github.com/networkx/networkx/blob/e941e68bfc521c7788461bb680041c09705a2b09/networkx/algorithms/simple_paths.py
# which is scheduled for release in networkx-2.5
# Until then the code is placed here

def all_simple_edge_paths(G, source, target, cutoff=None):
    """Generate lists of edges for all simple paths in the graph G from
    source to target.
    A simple path is a path with no repeated nodes.
    Parameters
    ----------
    G : NetworkX graph
    source : node
       Starting node for path
    target : nodes
       Single node or iterable of nodes at which to end path
    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.
    Returns
    -------
    path_generator: generator
       A generator that produces lists of simple paths.  If there are no paths
       between the source and target within the given cutoff the generator
       produces no output.
       For multigraphs, the list of edges have elements of the form `(u,v,k)`.
       Where `k` corresponds to the edge key.
    Examples
    --------
    Print the simple path edges of a Graph::
        >>> import networkx as nx
        >>>
        >>> g = nx.Graph()
        >>> g.add_edge('1', '2')
        >>> g.add_edge('2', '4')
        >>> g.add_edge('1', '3')
        >>> g.add_edge('3', '4')
        >>>
        >>> for path in sorted(nx.all_simple_edge_paths(g, '1', '4')):
        ...     print(path)
        [('1', '2'), ('2', '4')]
        [('1', '3'), ('3', '4')]
    Print the simple path edges of a MultiGraph. Returned edges come with
    their associated keys::
        >>> mg = nx.MultiGraph()
        >>> mg.add_edge(1,2,key='k0')
        'k0'
        >>> mg.add_edge(1,2,key='k1')
        'k1'
        >>> mg.add_edge(2,3,key='k0')
        'k0'
        >>>
        >>> for path in sorted(nx.all_simple_edge_paths(mg, 1, 3)):
        ...     print(path)
        [(1, 2, 'k0'), (2, 3, 'k0')]
        [(1, 2, 'k1'), (2, 3, 'k0')]
    Notes
    -----
    This algorithm uses a modified depth-first search to generate the
    paths [1]_.  A single path can be found in $O(V+E)$ time but the
    number of simple paths in a graph can be very large, e.g. $O(n!)$ in
    the complete graph of order $n$.
    References
    ----------
    .. [1] R. Sedgewick, "Algorithms in C, Part 5: Graph Algorithms",
       Addison Wesley Professional, 3rd ed., 2001.
    See Also
    --------
    all_shortest_paths, shortest_path, all_simple_paths
    """
    if source not in G:
        raise nx.NodeNotFound('source node %s not in graph' % source)
    if target in G:
        targets = {target}
    else:
        try:
            targets = set(target)
        except TypeError:
            raise nx.NodeNotFound('target node %s not in graph' % target)
    if source in targets:
        return []
    if cutoff is None:
        cutoff = len(G) - 1
    if cutoff < 1:
        return []
    if G.is_multigraph():
        for simp_path in _all_simple_edge_paths_multigraph(G, source, targets,
                                                           cutoff):
            yield simp_path
    else:
        for simp_path in _all_simple_paths_graph(G, source, targets, cutoff):
            yield list(zip(simp_path[:-1], simp_path[1:]))


def _all_simple_edge_paths_multigraph(G, source, targets, cutoff):
    if not cutoff or cutoff < 1:
        return []
    visited = [source]
    stack = [iter(G.edges(source, keys=True))]

    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.pop()
        elif len(visited) < cutoff:
            if child[1] in targets:
                yield visited[1:] + [child]
            elif child[1] not in [v[0] for v in visited[1:]]:
                visited.append(child)
                stack.append(iter(G.edges(child[1], keys=True)))
        else: #len(visited) == cutoff:
            for (u,v,k) in [child]+list(children):
                if v in targets:
                    yield visited[1:] + [(u,v,k)]

            stack.pop()
            visited.pop()
