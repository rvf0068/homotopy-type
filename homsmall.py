"""
Calculates the homotopy type of graphs and their clique graphs
"""
import argparse
import random
import re
import mogutda
import networkx as nx
from pycliques.simplicial import (
    SimplicialComplex,
    clique_complex)
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import (
    has_dominated_vertex,
    complete_s_collapse,
    complete_s_collapse_edges)
from pycliques.cliques import clique_graph as k
from pycliques.helly import is_clique_helly
from pycliques.lists import list_graphs
from pycliques.surfaces import open_neighborhood


def simplify_ht(graph):
    """Simplifies the graph for homotopy type purposes"""
    v_graph = complete_s_collapse(graph)
    ev_graph = complete_s_collapse_edges(v_graph)
    vev_graph = complete_s_collapse(ev_graph)
    return vev_graph


def _read_dong(dong):
    """Converts the set given by dong_matching into a TeX string"""
    n_critical = len(dong)
    if n_critical == 0:
        return (True, "Contractible")
    else:
        list_dong = list(dong)
        dimension = len(list_dong[0])
        if n_critical == 1:
            return (True, f"\\(S^{ {dimension-1} }\\)")
        else:
            for simp in list_dong:
                if len(simp) != dimension:
                    return (False, dong)
            return (True, f"\\(\\vee_{ {n_critical} }S^{ {dimension-1} }\\)")


def _read_betti_numbers(bettis):
    betti = "\\("
    for index, betti_number in enumerate(bettis):
        if betti_number != 0:
            if betti_number == 1:
                betti = betti + f"S^{ {index} }"
            else:
                betti = betti + f"\\vee_{ {betti_number} }S^{ {index} }"
            betti = betti + "\\vee "
    return betti[:-5] + "\\)"


def _shuff(lista):
    return random.sample(list(lista), len(lista))


def homotopy_type(graph):
    """Attempts to get a homotopy type using Dong's matching and vertex
    decomposability"""
    c_complex = clique_complex(graph)
    dong1 = c_complex.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    s_ht = simplify_ht(graph)
    c_complex2 = clique_complex(s_ht)
    dong2 = c_complex2.dong_matching()
    if _read_dong(dong2)[0]:
        return _read_dong(dong2)[1]
    attempts = 0
    while attempts < 5:
        attempts = attempts + 1
        dong3 = c_complex2.dong_matching(order_function=_shuff)
        if _read_dong(dong3)[0]:
            return _read_dong(dong3)[1]
    c_c = collapse(c_complex)
    if is_vertex_decomposable(c_c):
        return _read_betti_numbers(betti_numbers(s_ht))
    return betti_numbers(s_ht)


def dimension_sc(s_complex):
    """Returns the dimension of a simplicial complex"""
    dims = [len(f) for f in s_complex.face_set]
    dims.sort()
    return dims[-1] - 1


def betti_numbers(graph):
    """Computes the betti numbers of the complex of completes of a graph"""
    def simplify_list(bettis):
        simplified = bettis
        while len(simplified) > 0 and simplified[-1] == 0:
            del simplified[-1]
        return simplified
    the_simplices = [tuple(c) for c in nx.find_cliques(graph)]
    the_complex = mogutda.SimplicialComplex(simplices=the_simplices)
    dim = dimension_sc(the_complex)
    numbers = [the_complex.betti_number(i) for i in range(dim+1)]
    numbers[0] = numbers[0] - 1  # reduced betti number
    return simplify_list(numbers)


def max_degree(graph):
    """Returns the maximum degree of a graph"""
    degrees = dict(graph.degree())
    return max(degrees.values())


def is_shedding_vertex(the_complex, vertex):
    """Returns True if vertex is a shedding vertex of the complex"""
    link = the_complex.link(vertex)
    deletion = the_complex.deletion(vertex)
    facets_deletion = deletion.facet_set
    faces_link = link.facet_set
    intersection = facets_deletion & faces_link
    return len(intersection) == 0


def is_vertex_decomposable(the_complex):
    """Returns True if the complex is vertex decomposable"""
    if len(the_complex.facet_set) == 1:
        return True
    else:
        for vertex in the_complex.vertex_set:
            if (is_shedding_vertex(the_complex, vertex) and
                is_vertex_decomposable(the_complex.link(vertex)) and
                    is_vertex_decomposable(the_complex.deletion(vertex))):
                return True
        return False


def _facets_containing_simplex(simplicial_complex, simplex):
    facets = simplicial_complex.facet_set
    return len([facet for facet in facets if simplex.issubset(facet)])


def is_free_face(simplicial_complex, simplex):
    """Returns True if simplex is a free face of simplicial_complex"""
    return (not(simplex in simplicial_complex.facet_set) and
            _facets_containing_simplex(simplicial_complex, simplex) == 1)


def remove_simplex(simplicial_complex, simplex):
    """Returns the simplicial complex obtained by removing a simplex"""
    facet_containing = [f for f in simplicial_complex.facet_set
                        if simplex.issubset(f)][0]
    other_facets = simplicial_complex.facet_set - {facet_containing}
    good_facets = other_facets
    for vertex in simplex:
        good = True
        for g_facet in other_facets:
            if (facet_containing-{vertex}).issubset(g_facet):
                good = False
                break
        if good:
            good_facets = good_facets.union({facet_containing-{vertex}})
    vertices = set.union(*(set(s) for s in good_facets))
    return SimplicialComplex(vertices, facet_set=good_facets)


def has_free_face(simplicial_complex):
    """Returns a free face of simplicial_complex if it exists, None otherwise"""
    for face in simplicial_complex.all_simplices():
        if is_free_face(simplicial_complex, face):
            return face
    return None


def collapse(simplicial_complex, verbose=False):
    """Collapses simplicial_complex"""
    s_c = simplicial_complex
    all_done = False
    while not all_done:
        if len(s_c.vertex_set) in {0, 1}:
            return s_c
        free_face = has_free_face(s_c)
        if free_face is None:
            all_done = True
        elif len(free_face) == 0:
            return SimplicialComplex({}, {})
        else:
            if verbose:
                print(f"Free face: {free_face}")
            s_c = remove_simplex(s_c, free_face)
    return s_c


def _is_special_cutpoint(graph, vertex):
    """Returns True if vertex is a special cutpoint of graph"""
    neigh = open_neighborhood(graph, vertex)
    if neigh.size() == 0:
        for n_vertex in neigh:
            if graph.degree[n_vertex] == 1:
                return False
        return True
    return False


def _find_special_cutpoint(graph):
    """Returns a special cutpoint of graph if it exists, None otherwise"""
    for vertex in graph:
        if _is_special_cutpoint(graph, vertex):
            return vertex
    return None


def _h_type_clique_graph_cutpoint(graph, vertex):
    """Returns the homotopy type of the clique graph of graph with a special cutpoint"""
    h_graph = graph.subgraph(set(graph.nodes())-{vertex})
    k_h_type = homotopy_type(k(h_graph))
    s_neigh = open_neighborhood(graph, vertex).order()
    if k_h_type == "Contractible":
        if s_neigh == 2:
            return "\\(S^{1}\\)"
        else:
            return f"\\(\\vee_{ {s_neigh-1} }S^{ {1} }\\)"
    else:
        k_h_type = k_h_type[2:]
        if "S^{1}" in k_h_type:
            pat = r"\_\{(\d+)\}S\^\{1\}"
            m = re.search(pat, k_h_type)
            if m is None:
                # K(G-v) contains only one copy of S^1
                return f"\\(\\vee_{ {s_neigh} }" + k_h_type
            inds = m.span(1)
            newcadena = k_h_type[:inds[0]]+str(int(k_h_type[inds[0]: inds[1]])+1)+k_h_type[inds[1]:]
            return "\\("+newcadena
        if s_neigh == 2:
            return "\\(S^{1}\\vee " + k_h_type
        return f"\\(\\vee_{ {s_neigh-1} }S^{ {1} }\\vee " + k_h_type


HEADING = ("| index | order | max d | Helly | K Helly | HT G | HT KG |\n"
           "|-------+-------+-------+-------+---------+------+-------|\n")


def main():
    """Main function"""
    if args.order == 7:
        all_graphs = nx.graph_atlas_g()
        results = "homotopy_types_up_to_7.org"

        def conditions(graph):
            return (graph.order() > 1
                    and nx.is_connected(graph)
                    and not has_dominated_vertex(graph)
                    and max_degree(graph) >= 5)
    else:
        all_graphs = list_graphs(args.order)
        results = f"homotopy_types_{args.order}_{args.start}.org"

        def conditions(graph):
            return (not has_dominated_vertex(graph)
                    and max_degree(graph) >= 5)
    i = args.start
    with open(results, 'a', encoding="utf8") as the_file:
        the_file.write(HEADING)
        for graph in all_graphs[(args.start+1):]:
            i = i+1
            print("\r", end='')
            print(f"Currently on graph {i}", end='', flush=True)
            if conditions(graph):
                p_g = p(graph)
                h_g = homotopy_type(p_g)
                k_g = k(p_g, 23)
                if k_g is not None:
                    pkg = nx.convert_node_labels_to_integers(p(k_g))
                    hkg = homotopy_type(pkg)
                    c_v = _find_special_cutpoint(graph)
                    if c_v is not None:
                        hkg = _h_type_clique_graph_cutpoint(p_g, c_v)
                    else:
                        hkg = homotopy_type(pkg)
                    is_helly = is_clique_helly(graph)
                    is_k_helly = is_clique_helly(pkg)
                    if (not (is_helly and
                             ("S^{1}" in h_g or "S^{1}" in hkg)
                             )) and not (
                                is_k_helly and h_g == hkg and "S^{1}" in h_g):
                        the_file.write("|" + str(i) +
                                       "|" + str(p_g.order()) +
                                       "|" + str(max_degree(graph)) +
                                       "|" + str(is_helly) +
                                       "|" + str(is_k_helly) +
                                       "|" + str(h_g) +
                                       "|" + str(hkg) +
                                       "|\n")
                else:
                    c_big = f"|{i}|Clique graph has at least 23 vertices|||||\n"
                    the_file.write(c_big)
    print("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("order", type=int, help="Order of graphs examined")
    parser.add_argument("--start", type=int, help="Index to start", default=-1)
    args = parser.parse_args()

    main()
