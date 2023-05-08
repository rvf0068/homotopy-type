"""
Checks the homotopy type of all connected graphs with 8 vertices
"""
import networkx as nx
from pycliques.simplicial import clique_complex
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import has_dominated_vertex
from pycliques.cliques import clique_graph as k
from pycliques.lists import list_graphs


all_graphs = list_graphs(8)
RESULTS = 'homotopy_8.org'


def main():
    """Main function"""
    i = -1
    with open(RESULTS, 'a', encoding="utf8") as the_file:
        the_file.write("| index | order | HT G | HT KG |\n")
        the_file.write("|-------+-------+------+-------|\n")
        for graph in all_graphs:
            i = i+1
            if graph.order() > 1 and nx.is_connected(graph) and not has_dominated_vertex(graph):
                pared_graph = p(graph)
                h_g = clique_complex(pared_graph).dong_matching()
                pkg = nx.convert_node_labels_to_integers(p(k(pared_graph)))
                hkg = clique_complex(pkg).dong_matching()
                the_file.write(f"|{i}|{pared_graph.order()}|{h_g}|{hkg}|\n")


if __name__ == '__main__':
    main()
