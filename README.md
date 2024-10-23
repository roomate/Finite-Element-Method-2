# Finite Element Method: Analysis and approximation of PDEs

This repository stores the code of an introductory course to a non-conformal Finite Elements Method, a similar but yet fundamentally different from the first repo 'Finite-Elements-Method' on the matter. The resolution of parabolic and elliptic PDEs  turns differently as a 'non-conformal' approach is employed, in opposition to the conformal methods knwon as 'Galerkin' methods. More precisely, the term non-conformal signifies that the numerical solution is not admissible to the abstract variational problem, which is usually the case. Let's give a quick illustration.

On the left figure below is the numerical solution to Laplace equation in a square domain using this method.
A simple way to understand that the solution is 'non-conformal' is the existence of discontinuities points at the interface of the finite elements on the right figure below. Indeed, it is well known that functions in Sobolev space have to be continuous throughout any interfaces, which is obviously not the case here.

<img src="img/uh1.PNG" alt="drawing" width="370"/> <img src="img/uh1_coup.PNG" alt="drawing" width="350"/>
