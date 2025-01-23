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

def rhoSOR(A, omega):
    """
    Computes the spectral radius of the iteration matrix for the SOR method.
    """
    D = np.diag(np.diag(A))  # Diagonal part of A
    L = np.tril(A, -1)  # Lower triangular part of A
    U = np.triu(A, 1)  # Upper triangular part of A
    M = np.linalg.inv(D) @ (L + U)  # Iteration matrix for SOR
    spectral_radius = max(abs(eigvals(M)))
    return spectral_radius

def sor_method(A, b, omega=1.0, tol=None, max_its=100, x0=None):
    n = len(A)
    if tol is None:
        tol = [1e-8, 1e-8]
    elif np.isscalar(tol):
        tol = [tol, tol * tol]
    if x0 is None:
        x0 = np.zeros(n)
    
    norm_b = norm(b)
    x = x0.copy()
    r = b - A @ x
    norm_r = norm(r)
    its = 0
    
    # Compute spectral radius of the SOR iteration matrix
    spectral_radius = rhoSOR(A, omega)
    print(f"Spectral radius of the iteration matrix: {spectral_radius}")

    while (norm_r > tol[0] * norm_b + tol[1]) and (its < max_its):
        for i in range(n):
            sigma = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x[i] = (1 - omega) * x[i] + omega * (b[i] - sigma) / A[i, i]
        r = b - A @ x  # Compute residual
        norm_r = norm(r)
        its += 1

    if norm_r > tol[0] * norm_b + tol[1]:
        print("Warning: Unable to achieve the requested tolerance.")

    return x, its, norm_r

def solve_poisson_with_sor(N, omega=1.5, tol=1e-12, max_its=1000):
    # Grid setup
    xmax, xmin = 1, 0
    ymax, ymin = 1, 0
    h = (xmax - xmin) / (N - 1)
    x = np.linspace(xmin, xmax, N)
    y = np.linspace(ymin, ymax, N)
    X, Y = np.meshgrid(x, y)

    # Initialize boundary conditions (zero everywhere)
    U = np.zeros((N, N))

    # Exact solution and RHS
    u_exact = exact_solution(X, Y)
    f = rhs_function(X, Y)

    # Interior points
    num_interior = N - 2

    e = np.ones(num_interior)
    T = diags([e, -2 * e, e], [-1, 0, 1], shape=(num_interior, num_interior))
    I = eye(num_interior)
    laplacian = (kron(I, T) + kron(T, I)) / h**2

    # Flatten the RHS and exclude boundary points
    f_flat = f[1:-1, 1:-1].flatten()

    # Solve using SOR method
    u_flat, its, norm_r = sor_method(laplacian.toarray(), f_flat, omega=omega, tol=(tol, tol), max_its=max_its)

    # Reshape the solution to 2D and include boundary points
    u = np.zeros((N, N))
    u[1:-1, 1:-1] = u_flat.reshape((num_interior, num_interior))

    return X, Y, u, u_exact, its, norm_r

def main():
    # Grid size
    N = 9
    omega = 1.58  # Relaxation parameter

    # Solve Poisson equation with SOR method
    X, Y, u, u_exact, its, norm_r = solve_poisson_with_sor(N, omega)

    # Compute error
    error = np.max(np.abs(u - u_exact)) / np.max(np.abs(u_exact))
    print(f"Relative error: {error:.2e}")
    print(f"Iterations: {its}, Final Residual Norm: {norm_r:.2e}")

    # Plot numerical solution
    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, u, 50, cmap='viridis')
    plt.colorbar(label='Numerical Solution')
    plt.title(f'Numerical Solution (SOR)\nIterations: {its}, Residual Norm: {norm_r:.2e}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    # Plot exact solution
    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, u_exact, 50, cmap='viridis')
    plt.colorbar(label='Exact Solution')
    plt.title('Exact Solution')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    # Compare in 3D
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, u, cmap='viridis', edgecolor='k', alpha=0.8)
    ax.set_title('Numerical Solution (3D)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x, y)')
    plt.show()

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, u_exact, cmap='viridis', edgecolor='k', alpha=0.8)
    ax.set_title('Exact Solution (3D)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x, y)')
    plt.show()

if __name__ == "__main__":
    main()
