"""
Checks the homotopy type of all graphs up to 7 vertices
"""
import argparse
import networkx as nx
from pycliques.simplicial import clique_complex
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import has_dominated_vertex, complete_s_collapse, complete_s_collapse_edges
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


def homotopy_type(graph):
    """Attempts to get a homotopy type using Dong's matching"""
    c_complex = clique_complex(graph)
    dong1 = c_complex.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    else:
        c_complex = clique_complex(simplify_ht(graph))
        dong2 = c_complex.dong_matching()
        if _read_dong(dong2)[0]:
            return _read_dong(dong2)[1]
        else:
            return (dong1, dong2)


def main():
    """Main function"""
    if args.order == 7:
        all_graphs = nx.graph_atlas_g()
        results = "homotopy_types_up_to_7.org"
        def conditions(graph):
            return graph.order() > 1 and nx.is_connected(graph) and not has_dominated_vertex(graph)
    else:
        all_graphs = list_graphs(args.order)
        results = f"homotopy_types_{args.order}.org"
        def conditions(graph):
            return not has_dominated_vertex(graph)

    i = -1
    with open(results, 'a', encoding="utf8") as the_file:
        the_file.write("| index | order | Helly | HT G | HT KG |\n")
        the_file.write("|-------+-------+-------+------+-------|\n")
        for graph in all_graphs:
            i = i+1
            if conditions(graph):
                p_graph = p(graph)
                h_g = homotopy_type(p_graph)
                pkg = nx.convert_node_labels_to_integers(p(k(p_graph)))
                hkg = homotopy_type(pkg)
                is_helly = is_clique_helly(graph)
                if not (is_helly and ("S^1" in h_g or "S^1" in hkg)):
                    the_file.write(f"|{i}|{p_graph.order()}|{is_helly}|{h_g}|{hkg}|\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("order", type=int, help="Order of graphs examined")
    args = parser.parse_args()

    main()
