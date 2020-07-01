cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                            ccc
ccc  purpose: Inverse computation of a matrix                  ccc
ccc                                                            ccc   
ccc  Input and output varibles:                                ccc
ccc          input         n  size of matrix                   ccc
ccc          input         a  input matrix                     ccc
ccc          output        y  output inverse matrix            ccc
ccc                                                            ccc
ccc  Use:                                                      ccc
ccc          ludcmp        LU docomposition                    ccc
ccc          lubksb        find inverse of matrix by columns   ccc
ccc                                                            ccc
ccc  Warning: After matinv(a,n,y) is called, the input matrix  ccc
ccc           'a' would be changed (used as working space)     ccc
ccc           unless it is a diagnal matix. So if 'a' is       ccc
ccc           not a diagnal matrix and will be used after      ccc
ccc           calling matinv.f, 'a' should be saved before     ccc
ccc           calling matinv.f.                                ccc
ccc                                                            ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
     
     
      subroutine matinv(a,n,y)
c
      real a(n,n),y(n,n),indx(n)
      y=0.
      do i=1,n
        y(i,i)=1.
      end do
      call ludcmp(a,n,n,indx,d)
      do j=1,n
        call lubksb(a,n,n,indx,y(1,j))
      end do
c
      return
      end
c
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=100,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
c
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

