# RichardPoFlow
**Richard**s equation on **po**rous media for underground **flow** problems.

This project aims to explore uncertainty quantification (UQ) for ground flow problems.
The ground flow is modeled based on the Richards equation (h-based) using finite
difference method with Picard iteration.

### The code provides:
- Basic Richards equaion solver for 1d/2d/3d problems.
	- Accept Dirichlet and Neumann boundary conditions.
	- Accept user specified non-linear term function
	- Accept random permeability field input
- Fast optimized (vectorized) code for 1d Richards equation
- POD based reduced order model for 1d Richards equation
- DEIM POD based reduced order model for 1d Richards equation
- log-normal random field generator  
- log-normal random field approximating (dimension reduced) generator using KL theory.

### The code will provides:
- Data-driven surrogates for UQ problems
- Hybrid surrogates for UQ problems

### The code Folder:
- RandomField: functions for generating random fields
- 3dSolver: functions for 3d Richards equation
- 1dSloverAndRom: functions for 1d Richards equation and the reduced order model.
- Archive: code recycle bin just for the author.
- FasterOperation: Faster Operation library.
- EasyCode: Intuitive codes for dummy like me
(These codes may be computational slow and clumsy. They are not meant for practical use.
	I, however, believe that going through the inelegant but intuitive code is one of the best way
	to learn for *'engineering-driven'* folks like myself. So I decide to keep them
	and hop it might be useful for someone someday.)

###The code usage:

Run DEMO files at root for demonstrations.
Most function also comes with a test
or DEMO file for the purpose of demonstrations.



**The code is still under development.**

<!---
Thus, the code is developed in a way to be easily understood and modified.
Efficiency is not the priority and further vectorization is required if code efficiency is highly demanded.


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
-->
