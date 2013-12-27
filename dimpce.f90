Module dimpce
      implicit none

      integer :: id_proc, num_proc

      integer,parameter::MAXNMCS=1000000 

      integer::casemode
      integer::filenum
      integer::fctindx
      integer::dyncyccnt
      integer::OS,evlfnc,dynamics
      integer::probtype(20) ! 1=fixed 2=variable 
      integer::ndimt 

      double precision fmean,fvar,fstd
      double precision, dimension(20)   :: xavg,xstd,fmeanprime,fvarprime,xavgt,xstdt
      double precision,dimension(20,20)::fmeandbleprime,fvardbleprime
      
      double precision,dimension(20)::dat

      logical::mainprog,lhsdyn
      integer::OUUflag
      double precision,allocatable,dimension(:,:,:)::rmsemat
      character*60 :: outfile
      integer::runnum
      integer :: randomflag

      double precision :: diffloc2,difflocmax

      integer::fcnt,fgcnt,fghcnt

    End Module dimpce
