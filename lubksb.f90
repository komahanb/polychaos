SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  REAL*8:: a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL*8:: sum
  ii=0
  do  i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do  j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  end do
  do  i=n,1,-1
     sum=b(i)
     do  j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do

  return
end SUBROUTINE lubksb


!           C  (C) Copr. 1986-92 Numerical Recipes Software #,V41.04'21v.
