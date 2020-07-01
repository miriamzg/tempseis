c------------------------------------------------------------------------
c
c	Subroutine forward_modelling - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c	Note_1: Calls are made to : theo
c
c       Note_2: rmodel(1:nlayer): thickness
c               rmodel(nlayer+1:2*nlayer): velocity at upper interface
c               rmodel(2*nlayer+1:3*nlayer): velocity at lower interface
c               rmodel(3*nlayer+1:4*nlayer): Vp/Vs ratio
c
c-------------------------------------------------------------------------
c
	program test_W

	
	real		Amax, Amin, phi, Lmax, Lmin, pi, DTOR, RTOD
	real		rc_lon, rc_lat, depth, dip, strike, rc(3)
	real*4      W(3,3)
	real		dip_rad, strike_rad, phi_rad
	character*10	junk
	real		a(3,3), y(3,3)








	pi = 4.*atan(1.)
	DTOR = 2.*pi/360 
	RTOD = 1./ DTOR

	Lmax = 100.
	Lmin = 10.
	phi = 10.

	strike = 262
	dip = 26

	dip_rad = dip * DTOR
	strike_rad = strike * DTOR
	phi_rad = phi * DTOR


	write(*,*) Lmax, Lmin, phi, dip, strike
		
	

c---------------------------------------------------------------------------------------
c		FFT of the point source solution
c---------------------------------------------------------------------------------------
	call calc_W(Lmax, Lmin, phi_rad, dip_rad, strike_rad, W)



	return
	end


c -----------------------------------------------------------
	subroutine calc_W(Lmax, Lmin, phi, dip, strike, W)

	real*4 	M(3,3), S(3,3), W(3,3), I(3,3), R1(3,3), R2(3,3), R3(3,3), U(3,3), Sinv(3,3)
	real    phi, dip, strike, a, Lmax, Lmin, teta

	write(*,*) "Lmax, Lmin, phi, dip, strike ", Lmax, Lmin, phi, dip, strike

	pi = 4.*atan(1.)


c   Matrix containing the eigenvalues
	M(1,1) = Lmax
	M(1,2) = 0.0
	M(1,3) = 0.0
	M(2,1) = 0.0
	M(2,2) = Lmin
	M(2,3) = 0.0
	M(3,1) = 0.0
	M(3,2) = 0.0
	M(3,3) = 0.0

c   Identity matrix
	I(1,1) = 1.0
	I(1,2) = 0.0
	I(1,3) = 0.0
	I(2,1) = 0.0
	I(2,2) = 1.0
	I(2,3) = 0.0
	I(3,1) = 0.0
	I(3,2) = 0.0
	I(3,3) = 1.0

	teta = (pi/2.) - strike

c   Rotation matrix strike
	R1(1,1) = cos(teta)
	R1(1,2) = -sin(teta)
	R1(1,3) = 0.0
	R1(2,1) = sin(teta)
	R1(2,2) = cos(teta)
	R1(2,3) = 0.0
	R1(3,1) = 0.0
	R1(3,2) = 0.0
	R1(3,3) = 1.0

c   Rotation matrix dip
	R2(1,1) = 1.0
	R2(1,2) = 0.0
	R2(1,3) = 0.0
	R2(2,1) = 0.0
	R2(2,2) = cos(-dip)
	R2(2,3) = -sin(-dip)
	R2(3,1) = 0.0
	R2(3,2) = sin(-dip)
	R2(3,3) = cos(-dip)

c   Rotation matrix phi	
	R3(1,1) = cos(phi)
	R3(1,2) = -sin(phi)
	R3(1,3) = 0.0
	R3(2,1) = sin(phi)
	R3(2,2) = cos(phi)
	R3(2,3) = 0.0
	R3(3,1) = 0.0
	R3(3,2) = 0.0
	R3(3,3) = 1.0


	U = matmul(I,R1)
	U = matmul(U,R2)
	U = matmul(U,R3)

	S = U
	call matinv(S,3,Sinv)

	W = matmul(S,M)
	W = matmul(W,Sinv)


	return
	end 



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
      real a(n,n),y(n,n)
      integer indx(n)
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
















