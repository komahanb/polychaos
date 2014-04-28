mpif90 -fopenmp -g  -Wall -fbacktrace -pg -m64  -c main.f90
mpif90 -fopenmp -g  -Wall -fbacktrace -pg -m64  -o pc dimpce.o main.o higher.o LUroutines.o collsub.o srsmsub.o ludcmp.o lubksb.o svdcmp.o svbksb.o pythag.o randomroutines.o exactoutputfile.o sampdist.o dynsampdist.o evalcostf.o threebarcost.o setuprhs.o setupmat.o tecplot.o montecarlo.o mypoly.o mpi.o optimize.o BFGSroutines.o CalcstuffBFGS.o nieder.o sobol.o faure.o hammersley.o halton.o rmsebound.o Eulersolve.a libmir.a tapenade.a -L/usr/local/lib -lgsl -lgslcblas -lm -Wl,-rpath=.
 ----------- pc created ----------- 
