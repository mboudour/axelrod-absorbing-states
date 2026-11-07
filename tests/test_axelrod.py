import os, sys
import networkx as nx
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from counting import brute, formula, tree_count, cycle_count
from axsim import prepare, simulate, diagnose

SMALL = [nx.complete_graph(3), nx.path_graph(4), nx.star_graph(3), nx.cycle_graph(4),
         nx.Graph([(0, 1), (1, 2), (0, 2), (2, 3)]), nx.complete_graph(4)]

def test_count_formula_matches_brute_force():
    for G in SMALL:
        for F, q in [(2, 2), (2, 3), (3, 2)]:
            assert brute(G, F, q) == formula(G, F, q)

def test_q2_count_is_graph_independent():
    for G in SMALL:
        n = G.number_of_nodes()
        for F in (2, 3):
            assert brute(G, F, 2) == 2 ** (n + F - 1)

def test_closed_forms():
    for n in (3, 4, 5):
        for F, q in [(2, 2), (2, 3), (3, 2)]:
            if (q**F)**n <= 2e6:
                assert brute(nx.path_graph(n), F, q) == tree_count(n, F, q)
                assert brute(nx.cycle_graph(n), F, q) == cycle_count(n, F, q)

def test_simulation_reaches_absorbing_state_with_theorem_properties():
    G = nx.karate_club_graph()
    P = prepare(G)
    for seed in range(5):
        s, t, ev, Q0 = simulate(P, 3, 6, seed)
        d = diagnose(P[0], s)              # asserts every edge has overlap 0 or F
        assert d['split'][1] == 0          # no triangle with exactly one boundary edge
        assert d['quotient_proper']        # every feature properly colours the domain quotient
        assert d['quotient_clique'] <= 6
