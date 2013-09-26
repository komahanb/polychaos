# configuration -- change this to install in different directories,
#                  or if "make" doesn't work.
# prefix = /usr/local

#########################
#       TARGET          #
#########################

TARGET= pcestimate.a

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
FC = ifort

FFLAGS  = -r8 -O4 -openmp  #-traceback #-ftrapuv -check uninit -traceback #  -g -fpe3 # -traceback -debug all

# -zero -fpe0  -CB  -O0  -g3 -debug extended -ftrapuv -check all #-parallel # -check
# -openmp #-check

CFLAGS = -fPIC -Wall -O0 -g -I$(GSL_prefix)/include
LFLAGS =  -L$(GSL_prefix)/lib -lgsl -lgslcblas -lm
LIBS = -ldl -lstdc++



#export:
#    ar rvs pcestimate.a *.o

SRCS = dimpce.o main.o higher.o LUroutines.o collsub.o srsmsub.o ludcmp.o lubksb.o svdcmp.o svbksb.o pythag.o randomroutines.o exactoutputfile.o sampdist.o dynsampdist.o evalcostf.o setuprhs.o setupmat.o tecplot.o montecarlo.o mypoly.o mpi.o

OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) 
	   ar rvs $@ $(OBJS)
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
