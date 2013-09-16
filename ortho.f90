  subroutine ortho(dist,dimpc,x,pl,dpl,ddpl)
    implicit none

    real*8:: x,pl,dpl,ddpl
    integer :: dist,dimpc
    
    if (dist.eq.1) then
       call LEGENDRE(dimpc,X(jj),PL(jj,:),DPL(jj,:),ddpl(jj,:))
    else if (dist.eq.2) then
       call hermite(dimpc,X(jj),PL(jj,:),DPL(jj,:),ddpl(jj,:),2)
    else
       stop'Not implemented yet'
    end if
   
  end subroutine ortho
