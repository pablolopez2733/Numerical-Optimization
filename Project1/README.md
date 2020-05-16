# Project 1
## Minimization of the Goldstein Price function

In this project we implemented several trust region algorithms for the Golstein Price function. We used six main functions:
* pCauchy.m ->This function calculates the Cauchy point given a function and a point x0.
* ApGrad and ApHessian which approximate the gradient vector and the Hessian matrix for a function f at a given point x0.
* pDogLeg.m -> this is an implementation of the DogLeg method for solving the Subproblem of the trust region given a point  x0, the gradient at x0 and the Hessian at x0.
* mRC1 and mRC2 -> These functions both calculate a minimum of the function. RC1 uses the Cauchy point for solving the SPTR and RC2 uses the DogLeg Method. 

### References:
* Numerical optimization / Jorge Nocedal, Stephen J. Wright.(2006)

