#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 18:07:12 2025

@author: tommasomelotti
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye
from scipy.sparse.linalg import spsolve

def exact_solution(x, y):
    return np.sin(np.pi * x)**2 * np.sin(np.pi * y)**2

def rhs_function(x, y):
    return 2 * np.pi**2 * (
        np.cos(2 * np.pi * x) * np.sin(np.pi * y)**2 +
        np.cos(2 * np.pi * y) * np.sin(np.pi * x)**2
    )

def solve_poisson_sparse(N):
    # Grid setup
    h = 1.0 / (N - 1)  # Step size
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)

    # Exact solution and RHS
    u_exact = exact_solution(X, Y)
    f = rhs_function(X, Y)

    # Interior points count
    num_interior = (N - 2)

    # Sparse Laplacian using Kronecker products
    e = np.ones(num_interior)
    T = diags([e, -2 * e, e], [-1, 0, 1], shape=(num_interior, num_interior))
    I = eye(num_interior)
    laplacian = (kron(I, T) + kron(T, I)) / h**2

    # Flatten the RHS, exclude boundary points
    f_flat = f[1:-1, 1:-1].flatten()

    # Solve the linear system
    u_flat = spsolve(laplacian, f_flat)

    # Reshape solution to 2D and add boundary conditions
    u = np.zeros((N, N))
    u[1:-1, 1:-1] = u_flat.reshape((num_interior, num_interior))

    return X, Y, u, u_exact

def compute_error(u, u_exact):
    return np.max(np.abs(u - u_exact)) / np.max(np.abs(u_exact))

def main():
    grid_sizes = [10, 20, 40, 80]
    errors = []
    hs = []

    for N in grid_sizes:
        X, Y, u, u_exact = solve_poisson_sparse(N)
        error = compute_error(u, u_exact)
        errors.append(error)
        hs.append(1.0 / (N - 1))

    # Convergence plot
    plt.figure(figsize=(8, 6))
    plt.loglog(hs, errors, marker='o', label='Numerical Error')
    plt.loglog(hs, [h**2 for h in hs], linestyle='--', label='Expected $O(h^2)$')
    plt.xlabel('Grid spacing (h)')
    plt.ylabel('Relative Error (max norm)')
    plt.title('Convergence of Finite Difference Solution')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

    # Plot solution for the largest grid size
    X, Y, u, u_exact = solve_poisson_sparse(grid_sizes[-1])
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, u, 50, cmap='viridis')
    plt.colorbar(label='Numerical Solution')
    plt.title('Solution to Poisson Equation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

if __name__ == "__main__":
    main()
