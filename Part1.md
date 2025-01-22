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

