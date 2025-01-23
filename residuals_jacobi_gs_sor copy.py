import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye
from numpy.linalg import norm, eigvals

def exact_solution(x, y):
    return np.sin(np.pi * x)**2 * np.sin(np.pi * y)**2

def rhs_function(x, y):
    return 2 * np.pi**2 * (
        np.cos(2 * np.pi * x) * np.sin(np.pi * y)**2 +
        np.cos(2 * np.pi * y) * np.sin(np.pi * x)**2
    )

def jacobi_method(A, b, tol=(1e-6, 1e-6), max_its=100, x0=None):
    n = len(A)
    if x0 is None:
        x0 = np.zeros(n)
    
    tol_rel, tol_abs = tol
    norm_b = np.linalg.norm(b)
    x = x0
    r = b - np.dot(A, x)
    norm_r = np.linalg.norm(r)
    its = 0
    
    # Track residuals
    residuals = [norm_r]

    D = np.diag(np.diag(A))
    G = np.eye(n) - np.linalg.solve(D, A)

    while norm_r > tol_rel * norm_b + tol_abs and its < max_its:
        its += 1
        x = np.dot(G, x) + np.linalg.solve(D, b)
        r = b - np.dot(A, x)
        norm_r = np.linalg.norm(r)
        residuals.append(norm_r)
    
    return x, residuals

def gauss_seidel_method(A, b, tol=(1e-6, 1e-6), max_its=100, x0=None):
    n = len(A)
    if x0 is None:
        x0 = np.zeros(n)

    tol_rel, tol_abs = tol
    norm_b = np.linalg.norm(b)
    x = x0
    r = b - np.dot(A, x)
    norm_r = np.linalg.norm(r)
    its = 0

    # Track residuals
    residuals = [norm_r]

    P = np.tril(A)
    N = -np.triu(A, 1)

    while norm_r > tol_rel * norm_b + tol_abs and its < max_its:
        its +=1
        x = np.linalg.solve(P, N @ x + b)
        r = b - np.dot(A, x)
        norm_r = np.linalg.norm(r)
        residuals.append(norm_r)
    
    return x, residuals

def sor_method(A, b, omega=1.5, tol=(1e-6, 1e-6), max_its=100, x0=None):
    n = len(A)
    if x0 is None:
        x0 = np.zeros(n)

    tol_rel, tol_abs = tol
    norm_b = np.linalg.norm(b)
    x = x0
    r = b - np.dot(A, x)
    norm_r = np.linalg.norm(r)
    its = 0

    # Track residuals
    residuals = [norm_r]

    D = np.diag(np.diag(A))
    L = np.tril(A, -1)
    U = np.triu(A, 1)

    while norm_r > tol_rel * norm_b + tol_abs and its < max_its:
        for i in range(n):
            sigma = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x[i] = (1 - omega) * x[i] + omega * (b[i] - sigma) / A[i, i]
        r = b - np.dot(A, x)
        norm_r = np.linalg.norm(r)
        residuals.append(norm_r)
        its+=1
    
    return x, residuals

def solve_poisson(N, method="jacobi", omega=1.5, tol=1e-12, max_its=1000):
    xmax, xmin = 1, 0
    ymax, ymin = 1, 0
    h = (xmax - xmin) / (N - 1)
    x = np.linspace(xmin, xmax, N)
    y = np.linspace(ymin, ymax, N)
    X, Y = np.meshgrid(x, y)

    u_exact = exact_solution(X, Y)
    f = rhs_function(X, Y)

    num_interior = N - 2
    e = np.ones(num_interior)
    T = diags([e, -2 * e, e], [-1, 0, 1], shape=(num_interior, num_interior))
    I = eye(num_interior)
    laplacian = (kron(I, T) + kron(T, I)) / h**2

    f_flat = f[1:-1, 1:-1].flatten()

    if method == "jacobi":
        u_flat, residuals = jacobi_method(laplacian.toarray(), f_flat, tol=(tol, tol), max_its=max_its)
    elif method == "gauss_seidel":
        u_flat, residuals = gauss_seidel_method(laplacian.toarray(), f_flat, tol=(tol, tol), max_its=max_its)
    elif method == "sor":
        u_flat, residuals = sor_method(laplacian.toarray(), f_flat, omega=omega, tol=(tol, tol), max_its=max_its)

    u = np.zeros((N, N))
    u[1:-1, 1:-1] = u_flat.reshape((num_interior, num_interior))

    return X, Y, u, u_exact, residuals

def plot_residuals(residuals_jacobi, residuals_gs, residuals_sor):
    plt.figure(figsize=(10, 6))
    plt.semilogy(residuals_jacobi, label="Jacobi")
    plt.semilogy(residuals_gs, label="Gauss-Seidel")
    plt.semilogy(residuals_sor, label="SOR")
    plt.xlabel('Iteration')
    plt.ylabel('Residual Norm')
    plt.legend()
    plt.title('Residual Norm vs Iteration for Jacobi, Gauss-Seidel, and SOR')
    plt.grid(True)
    plt.show()

def main():
    N = 9  # Grid size
    omega = 1.5  # Relaxation parameter for SOR

    # Solve Poisson equation with each method 
    _, _, _, _, residuals_jacobi = solve_poisson(N, method="jacobi")
    _, _, _, _, residuals_gs = solve_poisson(N, method="gauss_seidel")
    _, _, _, _, residuals_sor = solve_poisson(N, method="sor", omega=omega)

    # Plot residuals for comparison
    plot_residuals(residuals_jacobi, residuals_gs, residuals_sor)

if __name__ == "__main__":
    main()
