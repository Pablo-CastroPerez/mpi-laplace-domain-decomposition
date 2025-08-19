# MPI Laplace Solver Framework

This project provides a **reusable parallel (MPI) framework** to solve the 2D Laplace equation.  
It uses a **geometric domain decomposition pattern** to distribute the computational domain across processes, and includes a modular design that separates:

- **Solver kernels** (e.g. Jacobi, Jacobi-SOR, Gauss–Seidel, Gauss–Seidel-SOR)  
- **Halo exchange routines** for nearest-neighbour communication  
- **Iteration control and convergence monitoring**

The framework is designed to be extendable: new solvers or alternative halo exchange strategies can be integrated with minimal changes to the parallel infrastructure.  

## Code Structure

### `main.c`
Entry point of the program.  
Responsible for:
- Initialising and finalising MPI.  
- Parsing command-line arguments (global grid sizes, convergence accuracy, maximum iterations).  
- Determining the local domain size per MPI process.  
- Selecting and launching the desired solver via `run_solver`.  


### `solvers.c`
Implements the numerical solvers for the 2D Laplace equation.  
Provides both memory management (init/finalise functions) and the iterative update rules for different methods:

- **Jacobi Method** (`jacobi_solver`)  
  Classical fixed-point iteration: updates each grid point as the average of its neighbours.  

- **Jacobi with Successive Over-Relaxation (SOR)** (`jacobi_sor_solver`)  
  Adds a relaxation factor `W` to accelerate convergence (if chosen correctly).  

- **Gauss–Seidel Method** (`gauss_seidel_solver`)  
  Similar to Jacobi, but updates are performed in-place, using newly computed values immediately.  

- **Gauss–Seidel SOR** (`gauss_seidel_sor_solver`)  
  Combines Gauss–Seidel with over-relaxation for faster convergence.  
  

---

### `solvers.h`
Header file declaring the solver interface.  
It exposes:
- Initialisation and finalisation functions for Jacobi and Gauss–Seidel solvers.  
- Function prototypes for each solver variant (Jacobi, Jacobi-SOR, Gauss–Seidel, Gauss–Seidel-SOR).  

### `utils.c`
Implements the utility functions that glue together MPI parallelisation and the solver kernels.  
Contains:

- **`initialise`**  
  Sets up the grid with boundary conditions (`LEFT_VALUE` on the left, `RIGHT_VALUE` on the right) and zero elsewhere.  
  Ensures both solution arrays (`u_k`, `u_kp1`) start consistent.  

- **`run_solver`**  
  Orchestrates the entire iterative process:  
  1. Allocates memory and calls the chosen solver’s initialisation.  
  2. Computes the initial residual norm.  
  3. Iteratively calls the solver kernel, performs halo swaps, and checks convergence.  
  4. Reports progress every `REPORT_NORM_PERIOD` iterations.  
  5. Finalises the solver and reports runtime statistics.  

- **`perform_halo_swap`**  
  Exchanges boundary rows between neighbouring MPI ranks (ghost cell update).  
  Uses **non-blocking MPI send/receive** (`MPI_Isend`, `MPI_Irecv`) and `MPI_Waitall` for efficiency.  

- **`get_residual`**  
  Computes the residual norm of the current solution compared to the discrete Laplace operator.  
  Uses `MPI_Allreduce` to combine local contributions into the global norm.  

---

### `utils.h`
Header file declaring the utility functions.  
Contains:
- Grid initialisation (`initialise`)  
- Solver driver (`run_solver`)  
- Communication (`perform_halo_swap`)  
- Residual computation (`get_residual`)
- 
---

### `config.h`
Defines the **global configuration parameters** that control solver behaviour and boundary conditions.  
Includes:

- **Boundary values**  
  ```c
  #define LEFT_VALUE  1.0
  #define RIGHT_VALUE 10.0
  ```
  Fixed Dirichlet boundary conditions at the left and right sides of the domain.
  - **Iteration control**  
  ```c
  #define MAX_ITERATIONS 1000000
  #define REPORT_NORM_PERIOD 1000
  ```
  Maximum number of iterations and reporting frequency for the residual norm.
  - **Solver selection**  
  ```c
  #define SOLVER_TO_USE 
  ```
  Chooses which algorithm is run:
	•	0 → Jacobi
	•	1 → Jacobi with over-relaxation (SOR)
	•	2 → Gauss–Seidel
	•	3 → Gauss–Seidel with SOR

 **Relaxation parameter**  
  ```c
  #define W 1.4 
  ```
Used in the over-relaxed variants (jacobi_sor, gauss_seidel_sor).
Must satisfy 1 < W < 2 for convergence (but may diverge if chosen poorly).
  
---
## Building and Running

### Makefile
For compiling all source files with `mpicc`. 

### archer2.srun
An example Slurm script to submit the code to the compute nodes of Archer2:
```bash
  sbatch archer2.srun
  ```

  
