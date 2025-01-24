# Intro

# 1. Finite difference discretization

Meshgrid: \
![Untitled](https://github.com/user-attachments/assets/24520ca8-367d-4c71-bbb1-d263a8846ecb)

Exact solution compared with solution by dense matrix: \
![Untitled](https://github.com/user-attachments/assets/c2b30609-c7a0-4e2f-8924-da48d591001a)

Convergence plot: \
![Untitled](https://github.com/user-attachments/assets/4127f8c7-9e37-4b3c-9cbc-19f47f06736a)


Computational time: No problems for lower values but the higher the values the exponentially many issues arise

# 2. Solving the linear system
## 2.1 Direct methods

Solution of the poisson equation using a sparse matrix \
![Untitled-1](https://github.com/user-attachments/assets/483d6c52-9d65-41bd-82cc-077499cbe62b)

Exact solution: \
![Untitled](https://github.com/user-attachments/assets/c3bb18a4-f688-40eb-99bc-d16f652cd0cd)

Computational time of dense vs sparse methods: \
![Untitled](https://github.com/user-attachments/assets/5c775dc3-b5f7-4530-8c06-886fee44cabb)

## 2.2 Iterative methods 

Solution of the poisson equation using Jacobi \
![Untitled-1](https://github.com/user-attachments/assets/c693c1a6-a4b8-47e2-8b2e-2a77117e0676)
![Untitled](https://github.com/user-attachments/assets/b785bd7f-cbbc-46df-b945-9eaca749d958)

Solution of the poisson equation using Gauss-Seidel: \
![Untitled](https://github.com/user-attachments/assets/8a9282f8-ed42-40cd-9bf1-35e4cb9fc133)
![Untitled-1](https://github.com/user-attachments/assets/09d825ed-4a43-4178-8009-c2f7ced00557)

Solution of the poisson equation using SOR: \
![Untitled](https://github.com/user-attachments/assets/4d05a187-fb2f-4b13-958b-11b3cc61117d)
![Untitled-1](https://github.com/user-attachments/assets/f6e60703-ae8f-4f32-8a04-d569f63870b5)

Exact solution: \
![Untitled](https://github.com/user-attachments/assets/c3bb18a4-f688-40eb-99bc-d16f652cd0cd)
![Untitled-1](https://github.com/user-attachments/assets/3c982afa-408e-416d-8a4f-4938e86de54a)

Number of iterations vs Relaxation Parameter: \
![Untitled](https://github.com/user-attachments/assets/f70f58d2-2b95-4a64-ad74-3939fd32404b)

Residual Norm vs Iterations: \
![Untitled-1](https://github.com/user-attachments/assets/acbd64f4-24f9-46e3-ac7a-6d837283269c)

# 3. Extensions to the solver


![Untitled](https://github.com/user-attachments/assets/98660617-7795-410b-819d-ff5d5d0a103b) \
![Untitled-1](https://github.com/user-attachments/assets/b32e2aa9-7b7e-4ab8-be7d-f99ce9e9b235) \
![Untitled](https://github.com/user-attachments/assets/74d6026a-2ac3-444f-9e52-bc182f3b3ef0) \
![Untitled-1](https://github.com/user-attachments/assets/c463e7ad-3f5e-4382-b168-693a9cb93e21)




