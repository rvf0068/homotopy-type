"""
Calculates the homotopy type of graphs and their clique graphs
"""
import argparse
import random
import re
import timeit
import mogutda
import networkx as nx
from sympy import symbols, Poly, Add, Mul
from pycliques.simplicial import (
    Simplex,
    SimplicialComplex,
    clique_complex)
from pycliques.dominated import completely_pared_graph as p
from pycliques.dominated import (
    has_dominated_vertex,
    complete_s_collapse,
    complete_s_collapse_edges)
from pycliques.cliques import clique_graph as k
from pycliques.helly import is_clique_helly
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


def read_betti_numbers(bettis):
    if sum(bettis) == 0:
        return "Contractible"
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
    if graph.order() == 1:
        return "Contractible"
    graph = simplify_ht(graph)
    disconnected_complement = h_type_as_join_complement(graph)
    if disconnected_complement:
        return disconnected_complement
    star_c = h_type_using_star_cluster(graph)
    if star_c:
        return star_c
    spec_n = h_type_by_special_neigh(graph)
    if spec_n:
        return spec_n
    spec_e = h_type_by_special_edges(graph)
    if spec_e:
        return spec_e
    cut_p = h_type_by_cutpoints(graph)
    if cut_p:
        return cut_p
    c_complex = clique_complex(graph)
    dong1 = c_complex.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    s_ht = nx.convert_node_labels_to_integers(simplify_ht(graph))
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
    return homotopy_type_s_c(c_complex)


def special_vertex_in_s_c(s_c):
    for vertex in s_c.vertex_set:
        if s_c.link(vertex).dimension() == 0:
            return vertex
    return None


def h_type_s_c_by_special_vertex(s_c):
    s_c = collapse(s_c)
    vertex = special_vertex_in_s_c(s_c)
    if vertex is not None:
        s_c2 = s_c.deletion(vertex)
        h_type = homotopy_type_s_c(s_c2)
        s_neigh = len(s_c.link(vertex).vertex_set)
        return suspend_string_with_s1s(h_type, s_neigh)
    return False


def suspend_string_with_s1s(h_type, s_neigh):
    if h_type == "Contractible":
        if s_neigh == 2:
            return "\\(S^{1}\\)"
        else:
            return f"\\(\\vee_{ {s_neigh-1} }S^{ {1} }\\)"
    else:
        h_type = h_type[2:]
        if "S^{1}" in h_type:
            pat = r"\_\{(\d+)\}S\^\{1\}"
            m = re.search(pat, h_type)
            if m is None:
                # string contains only one copy of S^1
                return f"\\(\\vee_{ {s_neigh} }" + str(h_type)
            inds = m.span(1)
            newcadena = h_type[:inds[0]]+str(int(h_type[inds[0]: inds[1]])+s_neigh-1)+h_type[inds[1]:]
            return "\\("+newcadena
        if s_neigh == 2:
            return "\\(S^{1}\\vee " + str(h_type)
        return f"\\(\\vee_{ {s_neigh-1} }S^{ {1} }\\vee " + str(h_type)


def homotopy_type_s_c(s_c):
    c_c = collapse(s_c)
    dong1 = c_c.dong_matching()
    if _read_dong(dong1)[0]:
        return _read_dong(dong1)[1]
    attempts = 0
    while attempts < 5:
        attempts = attempts + 1
        dong3 = c_c.dong_matching(order_function=_shuff)
        if _read_dong(dong3)[0]:
            return _read_dong(dong3)[1]
    h_type = h_type_s_c_by_special_vertex(s_c)
    if h_type:
        return h_type
    if is_vertex_decomposable(c_c):
        return read_betti_numbers(betti_numbers_c(c_c))
    return betti_numbers_c(c_c)


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
            vertex = list(s_c.vertex_set)[0]
            return SimplicialComplex(vertex_set={vertex},
                                     facet_set={Simplex({vertex})})
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


def find_special_cutpoint(graph):
    """Returns a special cutpoint of graph if it exists, None otherwise"""
    for vertex in graph:
        if _is_special_cutpoint(graph, vertex):
            return vertex
    return None


def h_type_clique_graph_cutpoint(graph, vertex):
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
                return f"\\(\\vee_{ {s_neigh} }" + str(k_h_type)
            inds = m.span(1)
            newcadena = k_h_type[:inds[0]]+str(int(k_h_type[inds[0]: inds[1]])+s_neigh-1)+k_h_type[inds[1]:]
            return "\\("+newcadena
        if s_neigh == 2:
            return "\\(S^{1}\\vee " + str(k_h_type)
        return f"\\(\\vee_{ {s_neigh-1} }S^{ {1} }\\vee " + str(k_h_type)


def star(s_complex, vertex):
    new_facets = {f for f in s_complex.facet_set if vertex in f}
    vertices = set.union(*(set(s) for s in new_facets))
    return SimplicialComplex(vertices, facet_set=new_facets)


def star_cluster(s_complex, simplex):
    new_facets = set()
    for vertex in simplex:
        new_facets = new_facets.union({f for f in s_complex.facet_set
                                       if vertex in f})
    vertices = set.union(*(set(s) for s in new_facets))
    return SimplicialComplex(vertices, facet_set=new_facets)


def intersection_complex(s_complex1, s_complex2):
    vertices = {x for x in s_complex1.vertex_set if s_complex2.function({x})}

    def _new_function(s):
        return s_complex1.function(s) and s_complex2.function(s)
    return SimplicialComplex(vertices, function=_new_function)


def h_type_as_suspension(graph):
    def increase_numbers(line, n):
        pattern = r"S\^{(\d+)}"

        def repl(match):
            number = int(match.group(1))
            increased_number = number + n
            return f"S^{ {increased_number} }"
        updated_line = re.sub(pattern, repl, line)
        return updated_line
    c_graph = nx.complement(graph)
    comps = [c_graph.subgraph(c).copy() for c in nx.connected_components(c_graph)]
    compsK2 = [s for s in comps if s.order() == 2]
    if len(compsK2) == 0:
        return None
    others = [s for s in comps if s.order() != 2]
    subg = nx.subgraph(graph, set.union(*(set(s) for s in others)))
    return increase_numbers(homotopy_type(subg), len(compsK2))


def betti_numbers_c(simplicial_complex):
    """Computes the betti numbers of the complex of completes of a graph"""
    def simplify_list(bettis):
        simplified = bettis
        while len(simplified) > 0 and simplified[-1] == 0:
            del simplified[-1]
        return simplified
    the_simplices = [tuple(c) for c in simplicial_complex.facet_set]
    the_complex = mogutda.SimplicialComplex(simplices=the_simplices)
    dim = dimension_sc(the_complex)
    numbers = [the_complex.betti_number(i) for i in range(dim+1)]
    numbers[0] = numbers[0] - 1  # reduced betti number
    return simplify_list(numbers)


def h_type_using_star_cluster(graph):
    if graph.order() == 1:
        return "Contractible"
    graph = nx.convert_node_labels_to_integers(graph)
    c_graph = nx.complement(graph)
    verts = [i for i in c_graph.nodes() if open_neighborhood(c_graph, i).size() == 0]
    if len(verts) == 0:
        return False
    else:
        vertex = verts[0]
        IG = clique_complex(graph)
        ST = star(IG, vertex)
        SC = star_cluster(IG, c_graph[vertex])
        int_c = intersection_complex(ST, SC)
        csc = collapse(int_c)
        h_type = homotopy_type_s_c(csc)
        if h_type:
            return read_betti_numbers([0]+betti_numbers_c(csc))
        if is_vertex_decomposable(csc):
            return read_betti_numbers([0]+betti_numbers_c(csc))
        return False


def list_to_polynomial(lst):
    x = symbols('x')
    terms = [i * x**p for p, i in enumerate(lst)]
    polynomial = Poly(sum(terms), x)
    return polynomial


def polynomial_to_list(polynomial):
    coefficients = polynomial.all_coeffs()
    coefficients.reverse()
    return coefficients


def h_type_as_join_complement(graph):
    x = symbols('x')
    c_graph = nx.complement(nx.convert_node_labels_to_integers(graph))
    comps = [c_graph.subgraph(c).copy() for c in nx.connected_components(c_graph)]
    if len(comps) > 1:
        compls = [clique_complex(nx.complement(s)) for s in comps]
        if all([is_vertex_decomposable(collapse(c)) for c in compls]):
            pols = [list_to_polynomial(betti_numbers_c(s)) for s in compls]
            expr_pols = [poly.as_expr() for poly in pols]
            the_list = [0]*(len(comps)-1)+polynomial_to_list(Poly(Mul(*expr_pols), x))
            return read_betti_numbers(the_list)
        return False
    return False


def is_complete_graph(graph):
    n = graph.order()
    return graph.size() == n*(n-1)/2


def is_disjoint_union_of_two_completes(graph):
    comps = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
    return len(comps) == 2 and all([is_complete_graph(h) for h in comps])


def h_type_by_special_neigh(graph):
    neighs = [(i, open_neighborhood(graph, i)) for i in graph.nodes()]
    # twok2 = nx.disjoint_union(nx.complete_graph(2), nx.complete_graph(2))
    # filt = [v for (v, nei) in neighs if nx.is_isomorphic(nei, twok2)]
    filt = [v for (v, nei) in neighs if is_disjoint_union_of_two_completes(nei)]
    filt2 = [v for v in filt if nx.is_connected(graph.subgraph(set(graph.nodes()-{v})))]
    if len(filt2) > 0:
        v = filt2[0]
        h = graph.subgraph(set(graph.nodes())-{v})
        h_type = homotopy_type(nx.convert_node_labels_to_integers(h))
        if h_type == "Contractible":
            return "\\(S^{1}\\)"
        else:
            h_type = h_type[2:]
            if "S^{1}" in h_type:
                pat = r"\_\{(\d+)\}S\^\{1\}"
                m = re.search(pat, h_type)
                if m is None:
                    # G-v contains only one copy of S^1
                    return "\\(\\vee_{2}" + str(h_type)
                inds = m.span(1)
                newcadena = h_type[:inds[0]]+str(int(h_type[inds[0]: inds[1]])+1)+h_type[inds[1]:]
                return "\\("+newcadena
            return "\\(S^{1}\\vee " + str(h_type)
    return False


def is_special_edge(graph, edge):
    n1 = open_neighborhood(graph, edge[0])
    n2 = open_neighborhood(graph, edge[1])
    inter = set(n1).intersection(set(n2))
    # avoid "threads" between components
    return len(inter) == 0 and graph.degree(edge[0]) > 2 and  graph.degree(edge[1]) > 2


def find_special_edges(graph):
    edges = graph.edges()
    return [e for e in edges if is_special_edge(graph, e)]


def h_type_by_special_edges(graph):
    x = symbols('x')
    graph = simplify_ht(graph)
    if nx.is_connected(graph):
        sp_edges = find_special_edges(graph)
        if sp_edges:
            c_graph = graph.copy()
            c_graph.remove_edges_from(sp_edges)
            comps = [c_graph.subgraph(c).copy() for c in nx.connected_components(c_graph)]
            if len(comps) == 2:
                compls = [clique_complex(s) for s in comps]
                if all([is_vertex_decomposable(collapse(c)) for c in compls]):
                    pols = [list_to_polynomial(betti_numbers_c(s)) for s in compls]
                    pols.append(list_to_polynomial([0, len(sp_edges)-1]))
                    expr_pols = [poly.as_expr() for poly in pols]
                    the_list = polynomial_to_list(Poly(Add(*expr_pols), x))
                    return read_betti_numbers(the_list)
            return False
        return False
    return False


def is_cutpoint(graph, vertex):
    if nx.is_connected(graph):
        c_graph = graph.copy()
        c_graph.remove_node(vertex)
        if not nx.is_connected(c_graph):
            return True
    return False


def find_cutpoints(graph):
    return [v for v in graph.nodes if is_cutpoint(graph, v)]


def h_type_by_cutpoints(graph):
    x = symbols('x')
    graph = simplify_ht(graph)
    if nx.is_connected(graph):
        cutpoints = find_cutpoints(graph)
        if cutpoints:
            vertex = cutpoints[0]
            c_graph = graph.copy()
            c_graph.remove_node(vertex)
            comps = [c.union({vertex}) for c in nx.connected_components(c_graph)]
            comps = [graph.subgraph(c) for c in comps]
            compls = [clique_complex(s) for s in comps]
            if all([is_vertex_decomposable(collapse(c)) for c in compls]):
                pols = [list_to_polynomial(betti_numbers_c(s)) for s in compls]
                expr_pols = [poly.as_expr() for poly in pols]
                the_list = polynomial_to_list(Poly(Add(*expr_pols), x))
                return read_betti_numbers(the_list)
            return False
        return False
    return False


HEADING = ("| index | order | max d | Helly | K Helly | HT G | HT KG |\n"
           "|-------+-------+-------+-------+---------+------+-------|\n")


def main():
    """Main function"""
    parts = args.filename.split('/')
    results = f"homotopy_types_{parts[-1]}.org"

    def conditions(graph):
        return (not has_dominated_vertex(graph)
                and max_degree(graph) >= 5)
    i = -1
    with open(results, 'a', encoding="utf8") as the_file:
        the_file.write(HEADING)
        graphs = nx.read_graph6(args.filename)
        for graph in graphs:
            i = i+1
            print("\r", end='')
            print(f"Currently on graph {i}", end='', flush=True)
            if conditions(graph):
                start_time = timeit.default_timer()
                p_g = p(graph)
                h_g = homotopy_type(p_g)
                k_g = k(p_g)
                if k_g is not None:
                    pkg = nx.convert_node_labels_to_integers(p(k_g))
                    c_v = find_special_cutpoint(graph)
                    if c_v is not None:
                        hkg = h_type_clique_graph_cutpoint(p_g, c_v)
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
                end_time = timeit.default_timer()
                print(f" Graph {i} took {end_time-start_time}")
        print("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of the file")
    args = parser.parse_args()

    main()
