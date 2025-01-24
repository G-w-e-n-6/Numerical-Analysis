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

Since the matrix of the linear system is a tri-diagonal matrix, it's useful to use a sparse representation where only the entries different from zero are stored. This representation allows us to save memory and computational cost. Usually, Gaussian elimination of an $$n \times n$$ matrix costs $$O(n^3)$$, but for a sparse matrix, we have $$O(n^2)$$ for the 2D Laplacian."


