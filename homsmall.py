"""
Checks the homotopy type of all graphs up to 7 vertices
"""
import argparse
import random
import mogutda
import networkx as nx
from pycliques.simplicial import clique_complex
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import (
    has_dominated_vertex,
    complete_s_collapse,
    complete_s_collapse_edges
    )
from pycliques.cliques import clique_graph as k
from pycliques.helly import is_clique_helly
from pycliques.lists import list_graphs


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


def _shuff(lista):
    return random.sample(list(lista), len(lista))


def homotopy_type(graph):
    """Attempts to get a homotopy type using Dong's matching"""
    c_complex = clique_complex(graph)
    dong1 = c_complex.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    else:
        s_ht = simplify_ht(graph)
        c_complex = clique_complex(s_ht)
        dong2 = c_complex.dong_matching()
        if _read_dong(dong2)[0]:
            return _read_dong(dong2)[1]
        else:
            attempts = 0
            while attempts < 5:
                attempts = attempts + 1
                dong3 = c_complex.dong_matching(order_function=_shuff)
                if _read_dong(dong3)[0]:
                    return _read_dong(dong3)[1]
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
    degrees = dict(graph.degree())
    return max(degrees.values())


def main():
    """Main function"""
    if args.order == 7:
        all_graphs = nx.graph_atlas_g()
        results = "homotopy_types_up_to_7.org"

        def conditions(graph):
            return (graph.order() > 1
                    and nx.is_connected(graph)
                    and not has_dominated_vertex(graph)
                    )
    else:
        all_graphs = list_graphs(args.order)
        results = f"homotopy_types_{args.order}.org"

        def conditions(graph):
            return (not has_dominated_vertex(graph)
                    and max_degree(graph) >= 6) 
    i = -1
    with open(results, 'a', encoding="utf8") as the_file:
        the_file.write("| index | order | Helly | K Helly | HT G | HT KG |\n")
        the_file.write("|-------+-------+-------+---------+------+-------|\n")
        for graph in all_graphs:
            i = i+1
            if conditions(graph):
                p_g = p(graph)
                h_g = homotopy_type(p_g)
                k_g = k(p_g, 23)
                if k_g is not None:
                    pkg = nx.convert_node_labels_to_integers(p(k_g))
                    hkg = homotopy_type(pkg)
                    is_helly = is_clique_helly(graph)
                    is_k_helly = is_clique_helly(pkg)
                    if not (is_helly and ("S^{1}" in h_g or "S^{1}" in hkg)):
                        the_file.write("|" + str(i) +
                                       "|" + str(p_g.order()) +
                                       "|" + str(is_helly) +
                                       "|" + str(is_k_helly) +
                                       "|" + str(h_g) +
                                       "|" + str(hkg) +
                                       "|\n")
                else:
                    the_file.write("|Clique graph has at least 23 vertices||||||\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("order", type=int, help="Order of graphs examined")
    args = parser.parse_args()

    main()
