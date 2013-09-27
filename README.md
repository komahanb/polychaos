=======================================================
Non-Intrusive Polynomial Chaos (Stochastic Collocation)
=======================================================

Author: Komahan Boopathy, University of Dayton, OH
Email : komahan.cool@gmail.com

=======================================================

A library that implements Polynomial Chaos (Regression). The targeted application is uncertainty quantification and using it for robust optimization.

1) PROBLEM SETUP
----------------

Multi-dimensional response surface is built using the input data. The user can control the choice of basis function. Currently the version has 

a) Legendre and

b) Hermite

orthogonal bases. 

2) SOLUTION
------------

The linear system is solved using scdcmp or ludcmp routines from "Numercial Recipes" depending upon whether the system is overdetermined or not. Please refer to their license information if you would like to taylor the code for your purposes.

3) POST PROCESSING:
--------------------

i)  Montecarlo simulation is carried out on the built response surface to give the statistics
 
ii) RMSE comparison between the surrogate and exact function can be carried out by changing the flag

iii)Tecplot output is written for funtions that have two variables or less.
