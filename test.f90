program main

  implicit none

  include 'mpif.h'
  integer :: NMCin,fctindxin,Casemode
  integer :: dim
  parameter (dim=3)

  double precision :: xavgin(dim),xstdin(dim),fmeanout,fvarout,fmeanprimeout(dim),fvarprimeout(dim)
 
  call MPI_START

  xavgin(:)=2.0
  xstdin(:)=0.15d0
  NMCin=100000
  fctindxin=0
  print *, 'hello'
  call PCestimate(dim,xavgin,xstdin,fmeanout,fvarout,fmeanprimeout,fvarprimeout,NMCin,fctindxin)

  !call Krigingestimate(ndimin,xavgin,ndimint,xstdin,fmeanout,fvarout,fmeanprimeout,fvarprimeout,NMCin,fctindxin)

  print *, fmeanout,fvarout

  call stop_all
end program main
