"""Shared paths: library code in src/, outputs in results/ (tracked) and paper/ (not tracked)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'src'))
RES = os.path.join(ROOT, 'results')
FIG = os.path.join(ROOT, 'paper', 'figures')
TAB = os.path.join(ROOT, 'paper', 'tables')
for d in (RES, FIG, TAB):
    os.makedirs(d, exist_ok=True)
