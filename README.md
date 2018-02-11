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

The linear system is solved using svdcmp or ludcmp routines from "Numercial Recipes" depending upon whether the system is overdetermined or not. Please refer to their license information if you would like to tailor the code for your purposes.

3) POST PROCESSING:
--------------------

i)  Montecarlo simulation is carried out on the built response surface to give the statistics
 
ii) RMSE comparison between the surrogate and exact function can be carried out by changing the flag

iii)Tecplot output is written for funtions that have two variables or less.


Citation:
---------

 K. Boopathy and M.P. Rumpfkeil,~\href{http://arc.aiaa.org/doi/abs/10.2514/1.J053064}
 {``Unified Framework for Training Point Selection and Error Estimation for Surrogate Models''}, AIAA
  Journal, Vol. 53, No. 1, pp. 215--234, 2015, DOI: 10.2514/1.J053064.


@article{Komahan2013c,
		  author      = "K. Boopathy and M. P. Rumpfkeil",
  		  TITLE        = "{Unified Framework for Training Point Selection and Error Estimation for Surrogate Models}",
                  journal     = "AIAA Journal",
		  volume      = "53, No. 1",
                  pages = "215--234",
                  doi = "http://dx.doi.org/10.2514/1.J053064",
                  year = "2015"
		  }
