! plots npts needed for pc for diff dimensions and order

program main
  implicit none
  integer*8::order
  integer*8::dim
  integer*8::npts

  integer*8::i,j,k,solver,stat,nterms,OS
  
  open(unit=43,file='nterms.dat')
  write(43,*) 'variables= "dim","order","f"'
  write(43,*) 'zone i=',10,' j=',10,' datapacking=point'
  solver=1
  OS=1
  do dim=1,10,1
     do stat=0,0 !F,FG,FGH
        do order=1,10 !order
           call combination(DIM+order,DIM,nterms) !get number of terms
           call getnpts(solver,stat,dim,nterms,OS,npts) !get number of points based on stat,solver,oversamp ratios
           write(43,*) dim, order, nterms
        end do
     end do
  end do
  close(43)

end program main



!+++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine getnpts(solver,stat,dim,nterms,oversamp,npts)

    implicit none

    integer*8,intent(in)  ::solver,stat,nterms,oversamp,dim
    integer*8,intent(out) ::npts


    if (stat.eq.0) then

       if (solver.eq.1) then !LU

          npts=nterms

       else !SV

          npts=oversamp*nterms

       end if

    else if (stat.eq.1) then

       if (solver.eq.1) then !LU
          npts=ceiling(real(nterms)/real(1+DIM))
       else

          npts=ceiling(real(oversamp*nterms)/real(1+DIM))
!          npts=npts+1

       end if

!       if (npts.lt.dim) npts=dim+1

    else if (stat.eq.2) then

       if (solver.eq.1) then !LU

          npts=ceiling(real(nterms)/(real(DIM)*real(1+DIM)/real(2.0d0)))

       else
          npts=ceiling(real(oversamp*nterms)/real(1+DIM+DIM*DIM))
!          npts=npts+dim
       end if

!       if (dim.eq.5) npts=npts+dim

!       if (npts.lt.dim) npts=dim+1

    else

       Stop'Unsupported Stat'

    end if


    return
  end subroutine getnpts


      subroutine combination(n,m,l)
      implicit none
      ! l = nCm = (n!)/(m!)/((n-m)!)
      integer*8, intent(in)  :: n,m
      integer*8, intent(out) :: l
      integer*8 :: i,i1,i2,i3,nbig,npet 

      if(n.le.0.or.m.le.0) stop'n<1 or m<1 in combination'
      if(n.lt.m) stop'n<m in combination'
      

      if(m.ne.n-m)then
         nbig = max(m,n-m)
         npet = min(m,n-m)
      else 
         nbig = m
         npet = n-m
      end if

      i1 = 1
      i2 = 1
      do i=nbig+1,n
         i1 = i1 * i
      end do
      do i=1,npet 
         i2 = i2 * i
      end do
      l = int( dble(i1)/dble(i2) )
  
      end subroutine combination
