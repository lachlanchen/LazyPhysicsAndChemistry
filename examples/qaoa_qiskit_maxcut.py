"""
QAOA Max-Cut with Qiskit (no Aer dependency).

Solves Max-Cut on a small 4-node ring graph using a p-layer QAOA.
Evaluation uses qiskit.quantum_info.Statevector to avoid external backends.
Also generates a heatmap figure (p=1) of expected cut vs (gamma, beta).

Run:
  python examples/qaoa_qiskit_maxcut.py
"""
from __future__ import annotations

import math
from typing import Iterable, List, Sequence, Tuple
import os

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


Edge = Tuple[int, int]


def maxcut_cost_from_bitstring(bits: Sequence[int], edges: Iterable[Edge]) -> float:
    # Cost is number of edges cut (endpoints differ)
    return sum(1.0 if bits[u] != bits[v] else 0.0 for (u, v) in edges)


def qaoa_circuit(n: int, edges: List[Edge], gammas: Sequence[float], betas: Sequence[float]) -> QuantumCircuit:
    p = len(gammas)
    assert p == len(betas)
    qc = QuantumCircuit(n)
    # Start in uniform superposition
    for i in range(n):
        qc.h(i)
    # Layers
    for layer in range(p):
        g = gammas[layer]
        b = betas[layer]
        # Cost unitary: exp(-i g sum_{(u,v) in E} Z_u Z_v)
        for (u, v) in edges:
            qc.rzz(2.0 * g, u, v)
        # Mixer unitary: exp(-i b sum_i X_i)
        for i in range(n):
            qc.rx(2.0 * b, i)
    return qc


def expected_cut_value_from_statevector(state: Statevector, edges: List[Edge], n: int) -> float:
    # Sum over computational basis probabilities weighted by cut size
    probs = state.probabilities()
    exp_cost = 0.0
    for idx, prob in enumerate(probs):
        if prob == 0.0:
            continue
        bits = [(idx >> i) & 1 for i in range(n)]
        exp_cost += prob * maxcut_cost_from_bitstring(bits, edges)
    return exp_cost


def objective(gammas: Sequence[float], betas: Sequence[float], n: int, edges: List[Edge]) -> float:
    qc = qaoa_circuit(n, edges, gammas, betas)
    sv = Statevector.from_instruction(qc)
    # We maximize cut value -> minimize negative expectation
    return -expected_cut_value_from_statevector(sv, edges, n)


def optimize_qaoa(n: int, edges: List[Edge], p: int = 1) -> Tuple[List[float], List[float], float]:
    # Try SciPy if available, else grid search fallback
    try:
        import numpy as np
        from scipy.optimize import minimize

        x0 = np.random.default_rng(42).uniform(0, math.pi, size=2 * p)

        def fun(x):
            gammas = x[:p]
            betas = x[p:]
            return objective(gammas, betas, n, edges)

        res = minimize(fun, x0, method="COBYLA", options={"maxiter": 200, "rhobeg": 0.5})
        xopt = res.x
        gammas = xopt[:p].tolist()
        betas = xopt[p:].tolist()
        val = -objective(gammas, betas, n, edges)
        return gammas, betas, val
    except Exception:
        # Simple coarse grid search for p=1
        best = (0.0, 0.0, -1.0)
        grid = [i * (math.pi / 16.0) for i in range(17)]
        for g in grid:
            for b in grid:
                val = -objective([g], [b], n, edges)
                if val > best[2]:
                    best = (g, b, val)
        return [best[0]], [best[1]], best[2]


def main():
    # 4-node ring graph
    n = 4
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    p = 1
    gammas, betas, exp_cut = optimize_qaoa(n, edges, p=p)
    max_possible = len(edges)
    print("QAOA (Qiskit) Max-Cut on 4-cycle")
    print(f"  p = {p}")
    print(f"  gammas = {gammas}")
    print(f"  betas  = {betas}")
    print(f"  expected cut value = {exp_cut:.4f} / {max_possible}")

    # Optional heatmap for p=1
    if p == 1 and plt is not None:
        import numpy as np
        ggrid = np.linspace(0, math.pi, 41)
        bgrid = np.linspace(0, math.pi, 41)
        Z = np.zeros((len(ggrid), len(bgrid)))
        for i, g in enumerate(ggrid):
            for j, b in enumerate(bgrid):
                qc = qaoa_circuit(n, edges, [g], [b])
                sv = Statevector.from_instruction(qc)
                Z[i, j] = expected_cut_value_from_statevector(sv, edges, n)

        # Locate max
        imax = np.unravel_index(np.argmax(Z), Z.shape)
        gopt = ggrid[imax[0]]
        bopt = bgrid[imax[1]]

        fig, ax = plt.subplots(figsize=(6, 5))
        hm = ax.imshow(Z, origin="lower", extent=[bgrid[0], bgrid[-1], ggrid[0], ggrid[-1]], aspect="auto")
        plt.colorbar(hm, ax=ax, label="Expected cut")
        ax.scatter([bopt], [gopt], c="w", s=30, marker="x", label="grid max")
        ax.set_xlabel("beta")
        ax.set_ylabel("gamma")
        ax.set_title("QAOA Max-Cut (4-cycle) expected cut")
        ax.legend(loc="upper right")

        out_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "figures"))
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "qaoa_qiskit_maxcut_heatmap.png")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        print(f"Saved heatmap: {out_path}")
    elif plt is None:
        print("matplotlib not available; skipping heatmap.")


if __name__ == "__main__":
    main()
