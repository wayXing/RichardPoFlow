# RichardPoFlow
Richards equation on porous media
using finite difference scheme.

This project aims to solve ground flow problems (based on Richard's equation with heterogeneous field input.)

One focus of this project is to develop general surrogate models (both data-driven and projection based reduced order model). The Richard's equation here also serves as a playground for these surrogate. Thus the code contains low efficient but easily readable version as well as more efficient version with less readability.

**The code is still under develop and not yet to be released.**


<!---
Thus, the code is developed in a way to be easily understood and modified.
Efficiency is not the priority and further vectorization is required if code efficiency is highly demanded.
-->

## Finished code.
- [x] permeability field (log-normal) generator based on KL decompositions.
- [x] Basic code for 1D/2D/3D domain (rectangular grid) with Dirichlet or  Neumann boundary conditions.
	* Proc: easy understanding and modifications
	* Cons: very slow! No vectorization/parallelization is yet developed for 2D/3D problems.
- [x] POD reduced order model for 1D problem.
- [x] hyper reduction for POD using discrete empirical interpolation. (for 1D problem).

## Ongoing code
- [ ] Data driven surrogate
- [ ] hybrid surrogate combing data-driven and projected based ROM.
