import networkx as nx
from pycliques.simplicial import clique_complex
from pycliques.dominated import pared_graph as p
from pycliques.cliques import clique_graph as k


all_graphs = nx.graph_atlas_g()
filelog = 'homotopy_types.org'
num_all_graphs = len(all_graphs)


def main():
    i = 0
    with open(filelog, 'a') as f:
        f.write("| index | order | HT G | HT KG |\n")
        f.write("|-------+-------+------+-------|\n")
        for g in all_graphs:
            pg = p(g)
            hg = clique_complex(pg).dong_matching()
            pkg = p(k(pg))
            hkg = clique_complex(pkg).dong_matching()
            f.write(f"|{i}|{pg.order()}|{hg}|{hkg}|")


if __name__ == '__main__':
    main()
