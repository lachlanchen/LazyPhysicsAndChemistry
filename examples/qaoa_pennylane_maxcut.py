"""
QAOA Max-Cut with PennyLane (no external chemistry deps).

Solves Max-Cut on a small 4-node ring graph using a p-layer QAOA.
Evaluation uses exact expectation on default.qubit (shots=None).

Run:
  python examples/qaoa_pennylane_maxcut.py
"""
from __future__ import annotations

import math
from typing import Iterable, List, Sequence, Tuple
import os

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

import pennylane as qml
from pennylane import numpy as np


Edge = Tuple[int, int]


def cost_hamiltonian(n: int, edges: Iterable[Edge]) -> qml.Hamiltonian:
    coeffs = []
    ops = []
    for (u, v) in edges:
        # 0.5 * (1 - Z_u Z_v)
        coeffs.append(0.5)
        ops.append(qml.Identity(0))  # I term (constant)
        coeffs.append(-0.5)
        ops.append(qml.PauliZ(u) @ qml.PauliZ(v))
    return qml.Hamiltonian(np.array(coeffs, dtype=float), ops)


def apply_cost_layer(edges: Iterable[Edge], gamma: float):
    for (u, v) in edges:
        qml.CNOT(wires=[u, v])
        qml.RZ(2.0 * gamma, wires=v)
        qml.CNOT(wires=[u, v])


def apply_mixer_layer(n: int, beta: float):
    for i in range(n):
        qml.RX(2.0 * beta, wires=i)


def build_qaoa_qnode(n: int, edges: List[Edge], p: int):
    H_cost = cost_hamiltonian(n, edges)
    dev = qml.device("default.qubit", wires=n, shots=None)

    @qml.qnode(dev)
    def circuit(params):
        # params shape: (p, 2) where [:,0]=gammas, [:,1]=betas
        gammas = params[:, 0]
        betas = params[:, 1]
        for i in range(n):
            qml.Hadamard(wires=i)
        for layer in range(p):
            apply_cost_layer(edges, gammas[layer])
            apply_mixer_layer(n, betas[layer])
        return qml.expval(H_cost)

    return circuit


def main():
    n = 4
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    p = 1
    circuit = build_qaoa_qnode(n, edges, p)
    params = np.random.default_rng(7).uniform(0, math.pi, size=(p, 2))
    opt = qml.GradientDescentOptimizer(stepsize=0.2)
    max_steps = 100
    history = []
    for step in range(max_steps):
        params, value = opt.step_and_cost(circuit, params)
        history.append(float(value))
        if (step + 1) % 20 == 0:
            print(f"Step {step+1:3d}  exp cut = {value:.6f}")
    final_val = circuit(params)
    print("QAOA (PennyLane) Max-Cut on 4-cycle")
    print(f"  p = {p}")
    print(f"  params = {params}")
    print(f"  expected cut value = {final_val:.6f} / {len(edges)}")

    if plt is not None and history:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(range(1, len(history) + 1), history, lw=2)
        ax.set_xlabel("Step")
        ax.set_ylabel("Expected cut")
        ax.set_title("PennyLane QAOA Max-Cut convergence")
        out_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "figures"))
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "qaoa_pennylane_maxcut_convergence.png")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        print(f"Saved convergence figure: {out_path}")
    elif plt is None:
        print("matplotlib not available; skipping convergence plot.")


if __name__ == "__main__":
    main()
