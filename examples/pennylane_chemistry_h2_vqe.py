"""
VQE for H2 (2-qubit reduced Hamiltonian) with PennyLane.

Uses a known 2-qubit Hamiltonian for H2 at ~0.735 Å (STO-3G, parity mapping
with two-qubit reduction). No external quantum chemistry packages required.

Reference Hamiltonian (from common tutorials, e.g., Qiskit VQE demos):
  H = c0*I + c1*Z0 + c2*Z1 + c3*Z0Z1 + c4*X0X1
with coefficients below.

Run:
  python examples/pennylane_chemistry_h2_vqe.py
"""
from __future__ import annotations

import os
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

import pennylane as qml
from pennylane import numpy as np


def h2_hamiltonian():
    # Coefficients for H2 at bond length ~0.735 Å
    c0 = -1.052373245772859
    c1 = 0.39793742484318045
    c2 = -0.39793742484318045
    c3 = -0.01128010425623538
    c4 = 0.18093119978423156

    coeffs = [c0, c1, c2, c3, c4]
    ops = [
        qml.Identity(0),
        qml.PauliZ(0),
        qml.PauliZ(1),
        qml.PauliZ(0) @ qml.PauliZ(1),
        qml.PauliX(0) @ qml.PauliX(1),
    ]
    return qml.Hamiltonian(np.array(coeffs, dtype=float), ops)


def ansatz_twolocal(params):
    # Simple 2-qubit hardware-efficient style ansatz
    # params shape: (2, 2) two rotation layers with entangling CNOTs in between
    qml.RY(params[0, 0], wires=0)
    qml.RY(params[0, 1], wires=1)
    qml.CNOT(wires=[0, 1])
    qml.RY(params[1, 0], wires=0)
    qml.RY(params[1, 1], wires=1)
    qml.CNOT(wires=[0, 1])


def main():
    H = h2_hamiltonian()
    dev = qml.device("default.qubit", wires=2, shots=None)

    @qml.qnode(dev)
    def energy_fn(params):
        ansatz_twolocal(params)
        return qml.expval(H)

    # Initialize and optimize
    params = np.array([[0.0, 0.0], [0.0, 0.0]], requires_grad=True)
    opt = qml.AdamOptimizer(stepsize=0.1)
    steps = 200
    history = []
    for i in range(steps):
        params, e = opt.step_and_cost(energy_fn, params)
        history.append(float(e))
        if (i + 1) % 20 == 0:
            print(f"Step {i+1:3d}  E = {e:.8f} Ha")

    e_final = energy_fn(params)
    print("PennyLane VQE for H2 (2-qubit reduced)")
    print(f"  Final energy: {e_final:.8f} Ha")
    print(f"  Parameters:\n{params}")

    if plt is not None and history:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(range(1, len(history) + 1), history, lw=2)
        ax.set_xlabel("Step")
        ax.set_ylabel("Energy (Ha)")
        ax.set_title("PennyLane VQE H2 convergence")
        out_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "figures"))
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "pennylane_h2_vqe_convergence.png")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        print(f"Saved convergence figure: {out_path}")
    elif plt is None:
        print("matplotlib not available; skipping convergence plot.")


if __name__ == "__main__":
    main()
