# Absorbing States of the Axelrod Model

Code for *Triadic and Assortative Structure of Absorbing States in the Axelrod Model* (M. Boudourides, 2026),
arXiv:2611.18719, [https://doi.org/10.48550/arXiv.2611.18719](https://doi.org/10.48550/arXiv.2611.18719).

**Interactive simulator.** https://mboudour.github.io/axelrod-absorbing-states/ (source: `docs/index.html`). It runs the
model on a chosen network and shows, side by side, the graph coloured by culture and the quotient graph of cultural domains,
with play, step and scrub controls, the edge and triangle composition over time, and the theorem checks at absorption.

**Setting.** Axelrod's model of cultural dissemination on a finite graph G = (V, E): each vertex carries a
culture in [q]^F, neighbours interact at a rate proportional to their overlap (number of shared features),
and an interaction copies one disagreeing feature (continuous-time, edge-driven version of Lanchier, 2012).
A configuration is absorbing when every edge is *interior* (overlap F) or *boundary* (overlap 0).

**What the code checks.**
1. The classification of absorbing configurations: a connected partition Π is a domain partition iff the
   quotient G/Π is properly q-colourable, and |A| = Σ_Π P(G/Π; q)^F (chromatic polynomials over connected
   partitions); closed forms for trees and cycles; |A| = 2^(n+F−1) for q = 2 on connected graphs.
2. Structure of absorbing states: no triangle has exactly one boundary edge; every feature properly
   colours the domain quotient; the assortativity chain r_f ≤ r_cult ≤ r_dom (Newman's categorical
   coefficient); open-wedge mass ω(G_I) ≤ ω(G) (two-path operator W = T + O of Boudourides,
   *J. Complex Netw.*, 2026; interior–boundary viewpoint of Boudourides, *DCDS-B*, 2026).
3. Dynamics: exact Gillespie simulation, validated against the generator, and the absorption-time bound
   E τ ≤ (F/2)(F n² − Q(σ₀)).

## Layout
```
src/axsim.py                       Gillespie simulator (numba) and diagnostics of absorbing states
src/counting.py                    brute-force enumeration, chromatic-sum formula, closed forms
src/workloads.py                   the nine graphs of Section 7 and the batch runner
scripts/01_benchmark.py            40 runs per graph and (F,q) ∈ {(2,2),(3,6),(3,12)} (1080 runs)
scripts/02_ws_sweep.py             Watts–Strogatz rewiring sweep (n=200, k=6, F=3, q=12)
scripts/03_enumeration.py          exhaustive enumeration on small graphs vs the counting formula
scripts/04_generator_validation.py first-transition law of the simulator vs exact rates (K2, P3, K3)
scripts/05_tables.py               LaTeX table bodies (written to paper/tables/, not tracked)
scripts/06_fig_karate.py           karate-club absorbing state and its domain quotient (pygraphviz)
scripts/07_fig_sweep.py            rewiring-sweep figure, cluster-bootstrap bands
scripts/08_lattice_animation.py    MP4 animation on a triangular lattice (not part of the paper)
docs/index.html                    interactive simulator (GitHub Pages)
tests/test_axelrod.py              unit tests (counting formula, q = 2 identity, closed forms, simulator)
results/                           JSON outputs used in the paper
```

## Reproduce
```bash
pip install -r requirements.txt pytest
cd scripts
python -m pytest -q ../tests
python 01_benchmark.py              # ~1 min
python 02_ws_sweep.py               # ~10 s
python 03_enumeration.py            # ~1 min
python 04_generator_validation.py   # ~2 min
python 05_tables.py
python 06_fig_karate.py
python 07_fig_sweep.py
python 08_lattice_animation.py      # optional, needs ffmpeg
```
Random seeds are fixed: in `01_benchmark.py` run j uses seed 1000+j; the random graphs are generated once
with fixed seeds; in `02_ws_sweep.py` graph realisation g uses seed 100+g. Rerunning the scripts reproduces
`results/*.json` exactly.

## Output (results/*.json)
`bench.json` and `sweep.json`: one record per run with the interior fraction `iota`, number of domains
`ndom` and cultures `ncul`, assortativity coefficients `r_f` (mean over features), `r_cul`, `r_dom`, the
triangle split `split` = numbers of triangles with 0, 1, 2, 3 boundary edges, open-wedge masses `omega_G`,
`omega_GI`, the colouring check `quotient_proper`, the clique number of the domain quotient
`quotient_clique`, absorption time `t`, and the total pairwise overlap `Q0`, `Qtau`.
`enum.json`: (graph, n, F, q, brute-force count, formula count). `genval.json`: total-variation distance
and mean waiting times of the simulator's first transition against the exact generator.

## License
MIT
