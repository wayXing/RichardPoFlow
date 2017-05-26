# RichardPoFlow
Richards equation on porous media
using finite difference

This project aims to solve ground flow problems (based on Richard's equation on heterogeneous field.)

One of the focus of this project is to develop general surrogate model (both by data driven and projection based reduced order model) to the
ground flow problems.  


<!---
Thus, the code is developed in a way to be easily understood and modified.
Efficiency is not the priority and further vectorization is required if code efficiency is highly demanded.
-->

## Finished code.
- [x]permeability fild (log-normal) generator based on KL decompositions.
- [x]Basic code for 1D/2D/3D domain (rectangular grid) with Dirichlet or  Neumann boundary conditions.
* Proc: easy understanding and modifications
* Cons: very slow! No vectorization/parallelization.

## Ongoing code
- [x] 1D case with POD
- [x] 1D case with POD-DEIM

## Further developments
- [ ] data driven surrogate
- [ ] hybrid surrogate wiht ML
