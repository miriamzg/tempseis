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
	program forward_modelling
c
        include 'rfi_param.inc'

	parameter     ( rad = 0.017453292 )
c	parameter    (maxdata = 5000)
c	parameter    (maxwave = 100)
	parameter    (nwave = 1)
c	parameter    (maxmoddim = 60)

c
	real*8		rmodel(maxmoddim)
	real*4          predicted_data_r(maxdata,maxwave),
     &                  predicted_data_i(maxdata,maxwave),
     &                  freq_array(maxdata, maxwave)
	integer		ndata(maxwave), nstations, nlines
	real		Amax, Amin, phi, Lmax, Lmin, pi, DTOR, RTOD
	character*200	command, filename
	character*10	stations(nwave)
	character*5	sta
	character*50	fault_file
	real		rc_lon, rc_lat, depth, dip, strike, rc(3)
	real		dip_rad, strike_rad, phi_rad
	character*10	junk
	real*4		dS_real(3,3,maxdata)
	real*4		dS_imag(3,3,maxdata)
	real*4		dS_freq(3,3,maxdata)
	real		a(3,3), y(3,3), W(3,3)
	real*4		freq(maxdata), point_source_real(maxdata), point_source_imag(maxdata)
	real*4		freq_E(maxdata), E_real(maxdata), E_imag(maxdata)
	real*4		Scalc_real(maxdata), Scalc_imag(maxdata)



c	fault_file = 'CMTSOLUTION_201604151625TEST_data.txt'
	fault_file = '201604151625A_TEST0_info.txt'
	open(12, file=fault_file, status="old")
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*)
	read(12,*) junk, rc_lon
	read(12,*) junk, rc_lat
	read(12,*) junk, rc_depth
	read(12,*) junk, dip
	read(12,*) junk, strike
	read(12,*) 
	read(12,*)
	read(12,*) junk, phi
	read(12,*) junk, Amax
	read(12,*) junk, Amin
	close(12)


	rc(1) = rc_lon
	rc(2) = rc_lat
	rc(3) = rc_dep

	Amax = rmodel(1)
c	Amin = rmodel(2)
c	phi = rmodel(3)

c	Amax = 20.
c	Amin = 5.

c	write(*,*) "Amax: ", Amax, " Amin: ", Amin, " phi: ", phi

	Lmax = (Amax**2)/4.
	Lmin = (Amin**2)/4.

	pi = 4.*atan(1.)
	DTOR = 2.*pi/360 
	RTOD = 1./ DTOR

	dip_rad = dip * DTOR
	strike_rad = strike * DTOR
	phi_rad = phi * DTOR

	do i=1,nwave
		sta = stations(i)
c		sta = "CO_0"
c----------------------------------------------
c		load kernels
c----------------------------------------------
		lw=lofw(sta)
		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dS2dx2',char(0)
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(1,1,j), dS_real(1,1,j), dS_imag(1,1,j)
		end do
		close(11)

		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dSdxdy',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(1,2,j), dS_real(1,2,j), dS_imag(1,2,j)
		end do
		close(11)

		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dSdxdz',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(1,3,j), dS_real(1,3,j), dS_imag(1,3,j)
		end do
		close(11)

		dS_freq(2,1,:) = dS_freq(1,2,:)
		dS_real(2,1,:) = dS_real(1,2,:)
		dS_imag(2,1,:) = dS_imag(1,2,:)
	
		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dS2dy2',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(2,2,j), dS_real(2,2,j), dS_imag(2,2,j)
		end do
		close(11)

		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dSdydz',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(2,3,j), dS_real(2,3,j), dS_imag(2,3,j)
		end do
		close(11)

		dS_freq(3,1,:) = dS_freq(1,3,:)
		dS_real(3,1,:) = dS_real(1,3,:)
		dS_imag(3,1,:) = dS_imag(1,3,:)

		dS_freq(3,2,:) = dS_freq(2,3,:)
		dS_real(3,2,:) = dS_real(2,3,:)
		dS_imag(3,2,:) = dS_imag(2,3,:)

		write(filename, '("kernels/",a,a9,a1)') sta(1:lw),'_Z_dS2dz2',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines
		read(11,*)
		do j=1,nlines
			read(11,*) dS_freq(3,3,j), dS_real(3,3,j), dS_imag(3,3,j)
		end do
		close(11)
c-----------------------------------------------------------------------------------
c		read point source solution
c-----------------------------------------------------------------------------------

		write(filename, '("point_source/",a,a6,a1)') sta(1:lw),'_Z_fft',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11,*) nlines
		ndata(i) = nlines
		read(11,*)
		do j=1, nlines
			read(11,*) freq(j), point_source_real(j), point_source_imag(j)
		end do
		close(11)

c---------------------------------------------------------------------------------------
c		FFT of the point source solution
c---------------------------------------------------------------------------------------
		call calc_W(Amax, Amin, phi_rad, dip_rad, strike_rad, W)
		
		do j=1,nlines
			E_real(j) = 0.0
			E_imag(j) = 0.0
			do k=1,3
				do l=1,3
					E_real(j) = E_real(j) +  W(k,l) * dS_real(k,l,j)
					E_imag(j) = E_imag(j) +  W(k,l) * dS_imag(k,l,j)
				end do
			end do
		end do



		do j=1,nlines
			predicted_data_r(j,i) = point_source_real(j) + 0.5 * E_real(j)
			predicted_data_i(j,i) = point_source_imag(j) + 0.5 * E_imag(j)
			freq_array(j,i) = freq(j)
		end do

		open(12, file=sta(1:lw)//"_predicted_real.txt", status="unknown")
		open(13, file=sta(1:lw)//"_predicted_imag.txt", status="unknown")
		do j=1,nlines
			write(12,*) freq_array(j,i), predicted_data_r(j,i)
			write(13,*) freq_array(j,i), predicted_data_i(j,i)
		end do
		close(12)
		close(13)
		


		
	end do	
	


	



	return
	end




c ----------------------------------------------------------------------------------------------------
c   CALC W matrix
c-----------------------------------------------------------------------------------------------------

	subroutine calc_W(Amax, Amin, phi, dip, strike, W)

	real 	M(3,3), R(3,3), W(3,3), I(3,3), R1(3,3), R2(3,3), R3(3,3), U(3,3), Rinv(3,3)
	real    phi, dip, strike, a, Amax, Amin, Lmax, Lmin, teta

	Lmax = (Amax/2.)**2
	Lmin = (Amin/2.)**2


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
	R = U

	W = matmul(R,M)
	call matinv(R,3,Rinv)
	W = matmul(W,Rinv)

	write(*,*) W(1,1), W(1,2), W(1,3)
	write(*,*) W(2,1), W(2,2), W(2,3)
	write(*,*) W(3,1), W(3,2), W(3,3)
	write(*,*) "--------------------"

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



c------------------------------------------
	integer function lofw(word)
c
	character*(*) word
c
	lw=len(word)
	k=0
	do i=1,lw
	  if (word(i:i).eq.' ') go to 99
	  k=k+1
	end do
99	lofw=k
c
	return
	end



