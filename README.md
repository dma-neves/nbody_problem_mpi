# Parallelizing the N-Body problem using MPI

![alt text](https://github.com/dma-neves/nbody_problem_mpi/blob/main/other/render_grid.jpg)

**Description**
  - In physics, the n-body problem is the problem of predicting the individual motions of a group of celestial objects interacting with each other gravitationally. Solving this problem has been motivated by the desire to understand the motions of the Sun, Moon, planets, and visible stars (Source: https://en.wikipedia.org/wiki/N-body_problem).
  - This project was made as an assignment for the *High Performance Computing* course, with the premise of studying the perfomance improvements gained from parallelizing a n-body algorithm using the MPI standard (implemented using OpenMPI). All of the algorithms and techniques used were based on those proposed by *P. Pacheco and M. Malensek. An Introduction to Parallel Programming. Elsevier Science, 2021*.
  - A complete explanation and performance analysis of these algorithms can be found in the full [report](https://github.com/dma-neves/nbody_problem_mpi/blob/main/report/report.pdf).

**Programs**
  - **n-body serial:** Serial/sequential aglorithm with both the basic and reduced solutions integrated (switched between the solutions using the `REDUCED` flag).
  - **n-body basic:** Parallel version of the basic solution. Every particle only calculates their receiving forces.
  - **n-body reduced:** Parallel version of the reduced solution. Every particle calculates receiving and symmetrical forces, reducing the number of iterations required.
  - **n-body render:** A particle trajectory renderer using the sequential reduced solution with some modifications (Not part of the assignment). A sequence of renderers of a 4 particle system can be seen above.

**Requirments**
  - Unix based OS.
  - gcc compiler.
  - Open MPI.
  - [CSFML library](https://www.sfml-dev.org/download/csfml/).
