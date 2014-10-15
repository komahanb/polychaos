# configuration -- change this to install in different directories,
#                  or if "make" doesn't work.
# prefix = /usr/local

#########################
#       TARGET          #
#########################

TARGET= pc

SUF90=f90
SUF77=f

.SUFFIXES: .f90 .f .o .mod

GSL_prefix = /usr/local

#########################
#      COMPILATION      #
#########################

FAD	= mpif90
F90	= mpif90
F77	= mpif77


CC = gcc
FC = gfortran

FFLAGS  = -fopenmp -g  -Wall -fbacktrace -pg -m64 #r8 -O4 -openmp # -fpe3 -parallel  #-traceback #-ftrapuv -check uninit -traceback #  -g -fpe3 # -traceback -debug all

# -zero -fpe0  -CB  -O0  -g3 -debug extended -ftrapuv -check all #-parallel # -check
# -openmp #-check

CFLAGS = -fPIC -Wall -O0 -g -I$(GSL_prefix)/include
LFLAGS =  -L$(GSL_prefix)/lib -lgsl -lgslcblas -lm
LIBS = -ldl -lstdc++



#export:
#    ar rvs pcestimate.a *.o

SRCS = dimpce.o dimcollsub.o main.o higher.o LUroutines.o ludcmp.o lubksb.o svdcmp.o svbksb.o pythag.o randomroutines.o exactoutputfile.o sampdist.o dynsampdist.o evalcostf.o threebarcost.o setuprhs.o setupmat.o tecplot.o montecarlo.o mypoly.o mpi.o optimize.o BFGSroutines.o CalcstuffBFGS.o nieder.o sobol.o faure.o hammersley.o halton.o rmsebound.o

OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) Eulersolve.a libmir.a tapenade.a
	$(F90) $(FFLAGS) -o $(TARGET) $(OBJS) Eulersolve.a libmir.a tapenade.a $(LFLAGS) -Wl,-rpath=.
	@echo " ----------- ${TARGET} created ----------- "

######################################
####### Compilation
######################################
%.o : %.mod

.$(SUF90).o:
	$(F90) $(FFLAGS) -c $<
.$(SUF77).o:
	$(F77) $(FFLAGS) -c $<

##################################
# Clean
##################################

clean:
	rm -f *.o *.so pc /output/output*
