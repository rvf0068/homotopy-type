"""
Checks the homotopy type of all graphs up to 7 vertices
"""
import networkx as nx
from pycliques.simplicial import clique_complex
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import has_dominated_vertex, complete_s_collapse, complete_s_collapse_edges
from pycliques.cliques import clique_graph as k


def simplify_ht(g):
    vg = complete_s_collapse(g)
    evg = complete_s_collapse_edges(vg)
    vevg = complete_s_collapse(evg)
    return vevg


def _read_dong(dong):
    n = len(dong)
    if n == 0:
        return (True, "Contractible")
    else:
        list_dong = list(dong)
        d = len(list_dong[0])
        if n == 1:
            return (True, f"\(S^{ {d-1} }\)")
        else:
            for simp in list_dong:
                if len(simp) != d:
                    return (False, dong)
            return (True, f"\(\\vee_{ {n} }S^{ {d-1} }\)")


def homotopy_type(g):
    cc = clique_complex(g)
    dong1 = cc.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    else:
        cc = clique_complex(simplify_ht(g))
        dong2 = cc.dong_matching()
        if _read_dong(dong2)[0]:
            return _read_dong(dong2)[1]
        else:
            return (dong1, dong2)


all_graphs = nx.graph_atlas_g()
RESULTS = 'homotopy_types.org'
num_all_graphs = len(all_graphs)


def main():
    """Main function"""
    i = -1
    with open(RESULTS, 'a', encoding="utf8") as the_file:
        the_file.write("| index | order | Helly | HT G | HT KG |\n")
        the_file.write("|-------+-------+-------+------+-------|\n")
        for graph in all_graphs:
            i = i+1
            if graph.order() > 1 and nx.is_connected(graph) and not has_dominated_vertex(graph):
                pared_graph = p(graph)
                h_g = homotopy_type(pared_graph)
                pkg = nx.convert_node_labels_to_integers(p(k(pared_graph)))
                hkg = homotopy_type(pkg)
                the_file.write(f"|{i}|{pared_graph.order()}|{is_clique_helly(g)}|{h_g}|{hkg}|\n")


if __name__ == '__main__':
    main()
