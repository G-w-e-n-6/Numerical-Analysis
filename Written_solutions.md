# 1.1 Designing the solver

2.
The finite difference formula for the second order derivative is as follows

$$\frac{\partial^2u}{\partial x^2} \approx \frac{u(x_{i+1}) - 2u(x_i) + u(x_{i-1})}{h^2}$$

h being the grid size.

Applying this formula to x and y directions results in the expressions:

$$\frac{\partial^2u}{\partial x^2} \approx \frac{u(x_{i+1},y_{j}) - 2u(x_i,y_j) + u(x_{i-1},y_j)}{h_x^2} ~~~ , ~~~ \frac{\partial^2u}{\partial y^2} \approx \frac{u(x_{i},y_{j+1}) - 2u(x_i,y_j) + u(x_i,y_{j-1})}{h_y^2}$$

$h_x$ and $h_y$ being the respective step size to each direction.

3.
The Laplacian formula is expressed as $\bigtriangleup = \frac{\partial^2u}{\partial x^2} + \frac{\partial^2u}{\partial y^2}$. So by adding the two partial derivatives and replacing $u(x_i,y_j)$ with $u_{i,j}$, the following result emerges:

$$ \bigtriangleup \approx \frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{h_x^2} + \frac{u_{i,j+1} - 2u_{i,j} + u_{i,j-1}}{h_y^2}$$

Since the grid for the $\Omega$ space is uniform, the step sizes in each direction are equivalent ($h_x = h_y$).

$$ \bigtriangleup \approx \frac{1}{h^2} (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j})$$


The equation is given as:

$$ D u_{j-1} + C u_j + D u_{j+1} = f_j, \quad 1 \leq j \leq n-1. $$

where

$$
C = 
\begin{pmatrix}
2 \left( \frac{1}{h_1^2} + \frac{1}{h_2^2} \right) & -\frac{1}{h_2^2} & 0 & \cdots & 0 \\
-\frac{1}{h_2^2} & 2 \left( \frac{1}{h_1^2} + \frac{1}{h_2^2} \right) & -\frac{1}{h_2^2} & \cdots & 0 \\
0 & -\frac{1}{h_2^2} & 2 \left( \frac{1}{h_1^2} + \frac{1}{h_2^2} \right) & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & -\frac{1}{h_2^2} \\
0 & 0 & 0 & -\frac{1}{h_2^2} & 2 \left( \frac{1}{h_1^2} + \frac{1}{h_2^2} \right)
\end{pmatrix},
$$

$$
D = 
\begin{pmatrix}
-\frac{1}{h_2^2} & 0 & 0 & \cdots & 0 \\
0 & -\frac{1}{h_2^2} & 0 & \cdots & 0 \\
0 & 0 & -\frac{1}{h_2^2} & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & 0 \\
0 & 0 & 0 & 0 & -\frac{1}{h_2^2}
\end{pmatrix},
$$

$$
f_j = f + b =
\begin{pmatrix}
f(x_1, y_j) + \frac{1}{h_1^2} \varphi(x_0, y_j) \\
f(x_2, y_j) \\
\vdots \\
f(x_{m-2}, y_j) \\
f(x_{m-1}, y_j) + \frac{1}{h_1^2} \varphi(x_m, y_j)
\end{pmatrix}.
$$

Equation can be further written as:

$$
\begin{pmatrix}
C & D & 0 & \cdots & 0 \\
D & C & D & \cdots & 0 \\
0 & D & C & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & D \\
0 & 0 & 0 & D & C
\end{pmatrix}
\begin{pmatrix}
u_1 \\
u_2 \\
\vdots \\
u_{n-2} \\
u_{n-1}
\end{pmatrix}=
\begin{pmatrix}
f_1 - D u_0 \\
f_2 \\
\vdots \\
f_{n-2} \\
f_{n-1} - D u_n
\end{pmatrix}.
$$

Since $u(x,1) = u(x, 0) = u(1, y) = u(0,y) = 0 $, the vector **b** is a zero-vector the size of $(N_y - 1)^2$. \
For a procedure with a non-zero boundary, the vector ressembles

$$ b =
\begin{pmatrix}
D u_0 \\
0 \\
\vdots \\
0 \\
D u_n
\end{pmatrix} $$

# 1.2 Validation of the implementation

The exact solution we are solving is given by:

$$
\
u(x, y) = \sin^2(\pi x) \sin^2(\pi y).
\
$$

The corresponding Poisson equation is:

$$
\
-\Delta u(x, y) = f(x, y),
\
$$ 


Substituting the exact solution \(u(x, y)\) into the Laplacian, the right-hand side \(f(x, y)\) is derived as:

$$
\
f(x, y) = 2\pi^2 \left( \cos(2\pi x) \sin^2(\pi y) + \cos(2\pi y) \sin^2(\pi x) \right).
\
$$

The exact solution is given by:

$$
u_{ex}(x, y) = \sin^2(\pi x) \sin^2(\pi y).
$$

We verify that $ u_{ex}(x, y) $ satisfies the boundary conditions:

1. At  x = 0 :

$$
u_{ex}(0, y) = \sin^2(\pi \cdot 0) \sin^2(\pi y) = 0.
$$

2. At  x = 1 :
   
$$
u_{ex}(1, y) = \sin^2(\pi \cdot 1) \sin^2(\pi y) = \sin^2(\pi) \sin^2(\pi y) = 0.
$$

4. At  y = 0 :

$$
u_{ex}(x, 0) = \sin^2(\pi x) \sin^2(\pi \cdot 0) = 0.
$$

5. At  y = 1 :

$$
u_{ex}(x, 1) = \sin^2(\pi x) \sin^2(\pi \cdot 1) = \sin^2(\pi x) \sin^2(\pi) = 0.
$$

Hence, $$u_{ex}(x, y)$$ is consistent with the boundary conditions:

$$
u(x, 0) = u(x, 1) = u(0, y) = u(1, y) = 0.
$$

# 2. Sovling the linear system

## 2.1 Direct methods

Since the matrix of the linear system is a tri-diagonal matrix, it's useful to use a sparse representation where only the entries different from zero are stored. This representation allows us to save memory and computational cost. Usually, Gaussian elimination of an $$n \times n$$ matrix costs $$O(n^3)$$, but for a sparse matrix, we have $$O(n^2)$$ for the 2D Laplacian.

## 2.2 Interative methods

The second derivative approximation on a uniform grid using finite differences can be expressed as:

$$
D_2 = \frac{1}{h^2}
\begin{bmatrix}
-2 & 1 & 0 & \cdots & 0 \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & 1 \\
0 & 0 & 0 & 1 & -2
\end{bmatrix},
$$

where $$h$$ is the grid spacing.

The eigenvalues $$\lambda_k$$ of the matrix $$D_2$$ are given by:

$$
\lambda_k = -\frac{4}{h^2} \sin^2\left(\frac{k \pi}{2n}\right),
$$

where $$k = 1, 2, \dots, n-1$$, and $$n$$ is the number of grid points.

The Chebyshev polynomials of the first kind, $$T_n(x)$$, satisfy the recurrence relation:

$$
T_0(x) = 1, \quad T_1(x) = x,
$$

$$
T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x), \quad n \geq 1.
$$

The eigenvalues of the second derivative matrix are related to the zeros of the Chebyshev polynomials $$T_n(x)$$ through the relation:

$$
x_k = \cos\left(\frac{(2k-1)\pi}{2n}\right), \quad k = 1, 2, \dots, n.
$$

In the 1D discrete case with Dirichlet boundary conditions, we are solving:

$$
\
\frac{v_{k+1} - 2v_k + v_{k-1}}{h^2} = \lambda v_k, \quad k = 1, \ldots, n, \quad v_0 = v_{n+1} = 0.
\
$$

Rearranging terms, we obtain:

$$
\
v_{k+1} = (2 + h^2 \lambda)v_k - v_{k-1}.
\
$$

Let $$2\alpha = 2 + h^2 \lambda$$. Assuming $$v_1 \neq 0$$, we can scale eigenvectors by any nonzero scalar, so scale $$v$$ such that $$v_1 = 1$$.

This leads to the recurrence relation:

$$
\
v_0 = 0, \quad v_1 = 1, \quad v_{k+1} = 2\alpha v_k - v_{k-1}.
\
$$

Considering \(\alpha\) as an indeterminate, we write:

$$
\
v_k = U_k(\alpha),
\
$$

where $$U_k$$ is the $$k$$-th Chebyshev polynomial of the second kind.

Since $$v_{n+1} = 0$$, we find:

$$
\
U_n(\alpha) = 0.
\
$$

Thus, the eigenvalues of the problem correspond to the zeros of the $$n$$-th Chebyshev polynomial of the second kind, with the relation:

$$
\
2\alpha = 2 + h^2 \lambda.
\
$$




The Kronecker product of two matrices affects their eigenvalues in the following way: if $$A$$ has eigenvalues $$\lambda_i$$ and $$B$$ has eigenvalues $$\mu_j$$, then the Kronecker product $$A \otimes B$$ has eigenvalues $$\lambda_i \mu_j$$.


Thus, for the 2D discrete Laplacian $$\Delta_h$$, the eigenvalues can be written as:

$$
\
\lambda_{i,j} = \lambda_i + \lambda_j,
\
$$

where $$\lambda_i$$ and $$\lambda_j$$ are the eigenvalues of the 1D discrete Laplacian. Therefore, for the 2D case, the eigenvalues $$\lambda_{i,j}$$ of $$\Delta_h$$ are:

$$
\
\lambda_{i,j} = 2 - 2 \cos\left(\frac{\pi i}{N + 1}\right) - 2 \cos\left(\frac{\pi j}{N + 1}\right),
\
$$

for $$i, j = 1, 2, \ldots, N$$.

This expression gives the eigenvalues of the 2D Laplacian in terms of the eigenvalues of the 1D Laplacian.

We expect these methods to converge because the Laplacian matrix is symmetric and semi-positive definite, since its eigenvalues are non-negative.

$$
\
\lambda_{i} = 2 - 2 \cos\left(\frac{\pi i}{N + 1}\right)
\
$$

The maximum of the absolute value of $$2 \cos\left(\frac{\pi i}{N + 1}\right)$$ is 2, so the eigenvalues are non-negative.

As the grid size increases, the smallest eigenvalue of the Laplacian approach zero.This increases the condition number of the matrix, leading to a spectral radius closer to 1, which slowes down convergence. 

The **convergence radius** of an iterative method is determined by the **spectral radius** of the iteration matrix \( T \), defined as the largest absolute value of its eigenvalues. 

For the iterative methods applied to the discrete Poisson problem, the convergence radius is:

$$
\
\rho(T) = \max |\lambda_i(T)|,
\
$$

where $$\lambda_i(T)$$ are the eigenvalues of the iteration matrix $$T$$.

###  Jacobi Method
The iteration matrix for Jacobi is:

$$
\
T_J = D^{-1}(L + U),
\
$$

where  $$D$$  is the diagonal part, and  $$L$$  and $$U$$  are the lower and upper triangular parts of the system matrix 

For the discrete Poisson problem, the convergence radius is:

$$
\
\rho(T_J) = \max \left| 1 - \frac{h^2 \lambda_i}{2} \right|,
\
$$


For larger grids, the smallest eigenvalues $$\lambda_i$$ of the Laplacian matrix approach zero, making $$\rho(T)$$ approach 1. This causes slower convergence for all methods, especially Jacobi and Gauss-Seidel.






### **Cost of an Iterative Solver**

The **cost of an iterative solver** depnd on two components:

1. **Cost per iteration**: The computational cost of performing a single iteration of the method.
2. **Number of iterations**: The total number of iterations required for the iterative solver to converge. 



#### Example: Jacobi Method
- The Jacobi method involves matrix-vector multiplications. The system matrix \( A \) (e.g., the discrete Laplacian) is sparse for PDEs like the Poisson equation, typically with \( O(N^2) \) unknowns in 2D for a grid of size $$N \times N$$.
- The sparsity of \( A \) means it has \( O(N^2) \) nonzero entries for a 5-point stencil.
- **Operations per iteration**:
  - A sparse matrix-vector multiplication costs \( O(N^2) \) operations.
  - Some additional vector operations (addition, scaling) are \( O(N^2) \).
- **Total cost per iteration**: \( O(N^2) \).


The **total cost** is the product of the cost per iteration and the number of iterations:

$$
\
\text{Total Cost} = \text{Cost per Iteration} \times \text{Number of Iterations}.
\
$$

For the Jacobi method:
- Cost per iteration: \( O(N^2) \).
- Number of iterations: Depends on $$\rho(T_J)$$ and the grid size $$N$$.

### **Optimal omega**

The optimal relaxation parameter $$\omega$$ for the Successive Over-Relaxation (SOR) method minimizes the spectral radius of the iteration matrix.


The convergence rate of SOR depends on the spectral radius $$\rho(T_{\text{SOR}})$$ where:

$$
\
T_{\text{SOR}} = (D - \omega L)^{-1} [(1 - \omega) D + \omega U],
\
$$

The smaller the spectral radius, the faster the convergence.



### **Optimal $$\omega_{\text{opt}}$$**
For the SOR method, it can be shown that:

$$
\
\omega_{\text{opt}} = \frac{2}{1 + \sin\left(\frac{\pi}{n+1}\right)}.
\
$$

# 3. Extensions to the solver
## 3.1  higher order finite difference methods
### 3.1.1  Richardson Extrapolation

The solution is denoted as $$u(h_1, h_2)$$.  

For the exact solution $$u(x, y)$$, we have:

**Theorem** Problem Definition:

$$
-\left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) = \frac{1}{12} \frac{\partial^4 u(x, y)}{\partial x^4}, \quad (x, y) \in \Omega,
$$

$$
v = 0, \quad (x, y) \in \Gamma.
$$

And:

$$
-\left( \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} \right) = \frac{1}{12} \frac{\partial^4 u(x, y)}{\partial y^4}, \quad (x, y) \in \Omega, 
$$

$$
w = 0, \quad (x, y) \in \Gamma. 
$$


**Existence of a Smooth Solution**

If a smooth solution exists, then we have:

$$
\max_{1 \leq i \leq m-1 \atop 1 \leq j \leq n-1} \left| u(x_i, y_j) - \left[ \frac{4}{3} u_{2i,2j} \left( \frac{h_1}{2}, \frac{h_2}{2} \right) - \frac{1}{3} u_{ij}(h_1, h_2) \right] \right| = O(h_1^4 + h_2^4),
$$

where $$h_1 = \frac{b-a}{m}$$, $$h_2 = \frac{d-c}{n}$$.

#### Proof

From above we can write:

$$
-(\delta_x^2 U_{ij} + \delta_y^2 U_{ij}) = f(x_i, y_j) - \frac{h_1^2}{12} \frac{\partial^4 u(x_i, y_j)}{\partial x^4} - \frac{h_2^2}{12} \frac{\partial^4 u(x_i, y_j)}{\partial y^4}
$$

$$
- \frac{h_1^4}{360} \frac{\partial^6 u(\xi_{ij}, y_j)}{\partial x^6} - \frac{h_2^4}{360} \frac{\partial^6 u(x_i, \eta_{ij})}{\partial y^6}, \quad (i,j) \in \omega,
$$

$$
u(x_i, y_j) = \varphi(x_i, y_j), \quad (i,j) \in \gamma,
$$

where $\xi_{ij} \in (x_i, x_{i+1})$$, $$\eta_{ij} \in (y_j, y_{j+1})$$.

The error equation can be rewritten as:

$$
-(\delta_x^2 e_{ij} + \delta_y^2 e_{ij}) = - \frac{h_1^2}{12} \frac{\partial^4 u(x_i, y_j)}{\partial x^4} - \frac{h_2^2}{12} \frac{\partial^4 u(x_i, y_j)}{\partial y^4}
$$

$$
- \frac{h_1^4}{360} \frac{\partial^6 u(\xi_{ij}, y_j)}{\partial x^6} - \frac{h_2^4}{360} \frac{\partial^6 u(x_i, \eta_{ij})}{\partial y^6}, \quad (i,j) \in \omega,

$$
e_{ij} = 0, \quad (i,j) \in \gamma.
$$

we construct the difference schemes:

$$
-(\delta_x^2 v_{ij} + \delta_y^2 v_{ij}) = \frac{1}{12} \frac{\partial^4 u(x_i, y_j)}{\partial x^4}, \quad (i,j) \in \omega,
\]

\[
v_{ij} = 0, \quad (i,j) \in \gamma,
\]

and

\[
-(\delta_x^2 w_{ij} + \delta_y^2 w_{ij}) = \frac{1}{12} \frac{\partial^4 u(x_i, y_j)}{\partial y^4}, \quad (i,j) \in \omega,
\]

\[
w_{ij} = 0, \quad (i,j) \in \gamma.
\]

Thus, by Theorem 2.5, we obtain:

\[
v(x_i, y_j) - v_{ij}(h_1, h_2) = O(h_1^2 + h_2^2), \quad (i,j) \in \omega,
\]

\[
w(x_i, y_j) - w_{ij}(h_1, h_2) = O(h_1^2 + h_2^2), \quad (i,j) \in \omega.
\]

