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
ca
	subroutine forward_modelling(
     &		rmodel, moddim, ndata, nwave, stations,  
     &		predicted_data, observed_data, 
     &		dS1, dS2, point_source, freq_obs, nlines, output_final)
c
        include 'rfi_param.inc'

c
	parameter     ( v60 = 8.043 )
	parameter     ( rad = 0.017453292 )
c
	integer		moddim
	real*8		rmodel(maxmoddim)
	real*4          predicted_data_r(maxdata,maxwave),
     &                  predicted_data_i(maxdata,maxwave)
	complex*8	predicted_data(maxdata,maxwave),
     &			observed_data(maxdata,maxwave)
	integer		ndata(maxwave), nstations
	real		Amax, Amin,  phi, Lmax, Lmin, pi, DTOR, RTOD, duration, Amin_Amax_ratio
	character*200	command, filename
	character*10	stations(maxwave)
	character*8	sta
	character*50	fault_file
	real		rc_lon, rc_lat, rc_depth, dip, strike, rc(3), vc(3), vel_abs
	real		rc_module, rc_angle, rc_angle_rad
	real		dip_rad, strike_rad, phi_rad, vel_angle, vel_angle_rad
	character*10	junk
	
	real		a(3,3), y(3,3), W(3,3)
	complex*8	W_imag(3,3)
	real*4		freq_obs(maxwave, maxdata)
	real*4		Scalc_real(maxdata), Scalc_imag(maxdata)
	complex*8	i_unit
	complex*8	AA(maxdata),BB(maxdata),CC(maxdata),DD(maxdata),EE(maxdata)
	complex*8	tensor_dot, dot, DD_dot
	real*4		tensor_dot_real, tensor_dot_imag
	logical		output_final, output_temp

	complex*8	dS1(maxwave, 3, maxdata)
	complex*8	dS2(maxwave, 3,3,maxdata)
	complex*8	point_source(maxwave, maxdata)
	integer		nlines(maxwave)

	i_unit = (0,1)


	output_temp = .false.



	n_param_lines = 6
c	rc_lon = rmodel(1)
c	rc_lat = rmodel(n_param_lines+1)
c	rc_depth = rmodel(2*n_param_lines+1)
	rc_module = rmodel(1)
	rc_angle = rmodel(n_param_lines+1)
	delta_t = rmodel(2)
	duration = rmodel(3)
	vel_abs = rmodel(4)
	vel_angle = rmodel(n_param_lines+4)
	Amax = rmodel(5)
	Amin_Amax_ratio = rmodel(n_param_lines+5)
	Amin = Amax * Amin_Amax_ratio
c---------------
c provvisorio
c	Amin = rmodel(n_param_lines+5)
c----------------

	phi = rmodel(2*n_param_lines+5)
	strike = rmodel(6)
	dip = rmodel(n_param_lines+6)



	pi = 4.*atan(1.)
	DTOR = 2.*pi/360 
	RTOD = 1./ DTOR

	rc_angle_rad = rc_angle * DTOR
	vel_angle_rad = vel_angle * DTOR
	dip_rad = dip * DTOR
	strike_rad = strike * DTOR
	phi_rad = phi * DTOR



c=============================
c	call v_given_par(Amax, Amin, vel_angle_rad, phi_rad, duration, vel_abs)
c	write(*,*) "V calc: ", vel_abs
c	stop
c=============================



	write(*,*) "-------------"
	write(*,*) "Spatial centroid module: ", rc_module, " angle: ", rc_angle
	write(*,*) "Temporal centroid: ",delta_t 
	write(*,*) "Duration: ", duration * sqrt(3.)
	write(*,*) "Velocity abs: ", vel_abs, " velocity angle: ", vel_angle
	write(*,*) "Amax: ", Amax * sqrt(3.), " Amin: ", Amin * sqrt(3.), " phi: ", phi
	write(*,*) "Strike: ", strike, " dip: ", dip 
	write(*,*) "-------------"


	write(123,*) "Spatial centroid module: ", rc_module, " angle: ", rc_angle
	write(123,*) "Temporal centroid: ",delta_t 
	write(123,*) "Duration: ", duration * sqrt(3.)
	write(123,*) "Velocity abs: ", vel_abs, " velocity angle: ", vel_angle
	write(123,*) "Amax: ", Amax * sqrt(3.), " Amin: ", Amin * sqrt(3.), " phi: ", phi
	write(123,*) "Strike: ", strike, " dip: ", dip 
	write(123,*) "***********************"








	call calc_vc(vel_abs, vel_angle_rad, strike_rad, dip_rad, vc)
	
	call calc_vc(rc_module, rc_angle_rad, strike_rad, dip_rad, rc)
	

c	rc(1) = rc_lon
c	rc(2) = rc_lat
c	rc(3) = rc_depth

c	write(*,*) delta_t

	do i=1,nwave
		sta = stations(i)
		lw=lofw(sta)
		
		call calc_W(Amax, Amin, phi_rad, dip_rad, strike_rad, W)

c		A = 1j * S1_freq * -tc  * S1_fft
c		C = - 0.5 * (S1_freq**2 / 2.)  * (Ttest/2.)**2 * S1_fft  
c		D = - 0.5 * 1j * S1_freq * (Ttest/2.)**2 * np.dot(vc, dS1)
c		E =   np.tensordot(W, dS2) / 2.

		do j=1,nlines(i)
			AA(j) = i_unit * freq_obs(i,j) * (-delta_t) * point_source(i,j) 
c--------------------------------------------------------------------------------------------
			BB(j) = 0.0
			do k=1,3
				BB(j) = BB(j) + rc(k) * dS1(i,k,j)
			end do
c--------------------------------------------------------------------------------------------
			CC(j) = -(freq_obs(i, j)**2 / 2.0) * (duration/2.0)**2 * point_source(i,j)
c--------------------------------------------------------------------------------------------
			DD_dot = 0.0
			do k=1,3
				DD_dot = DD_dot +  (vc(k)*dS1(i,k,j))
			end do
			DD(j) =  -i_unit * freq_obs(i,j) * (duration/2.0)**2 * DD_dot
c--------------------------------------------------------------------------------------------
			EE(j) = (0.0,0.0)
			tensor_dot = (0.0,0.0)
			do k=1,3
				do l=1,3
					tensor_dot = tensor_dot + W(k,l) * dS2(i,k,l,j)
				end do
			end do
			EE(j) =  0.5 * tensor_dot
c--------------------------------------------------------------------------------------------
c			write(*,*) "E ", j, point_source(j), tensor_dot, EE(j)
c			stop


			predicted_data(j,i) = point_source(i,j) + AA(j) + BB(j) + CC(j) + DD(j) + EE(j)
		end do	

c		open(666,file="predicted" , status="replace")
c		do j=1,nlines(i)
c			write(666, *) freq_obs(1,j), real(predicted_data(j,1)), imag(predicted_data(j,1))
c		end do
c		close(666)
c		stop

	end do




	if (output_temp) then
		do i=1,nwave
			sta = stations(i)
			lw=lofw(sta)
			write(filename, '("./rfi_files/TEMP/",a,a14,a1)') 
     &				sta(1:lw),'_predicted.asc',char(0)
			open(12, file=filename, status="replace")
			do j=1,nlines(i)
				write(12,*) freq_obs(i,j), real(predicted_data(j,i)), imag(predicted_data(j,i))
			end do	
			close(12)
			
		end do

		do i=1,nwave
			sta = stations(i)
			lw=lofw(sta)
			write(filename, '("./rfi_files/TEMP/",a,a13,a1)') 
     &				sta(1:lw),'_observed.asc',char(0)
			open(12, file=filename, status="replace")
			do j=1,nlines(i)
				write(12,*) freq_obs(i,j), real(observed_data(j,i)), imag(observed_data(j,i)) 
			end do	
			close(12)
			
		end do


	end if


	if (output_final) then
		do i=1,nwave
			sta = stations(i)
			lw=lofw(sta)
			write(filename, '("./rfi_files/FINAL/",a,a19,a1)') 
     &				sta(1:lw),'_best_predicted.asc',char(0)
			open(12, file=filename, status="replace")
			do j=1,nlines(i)
				write(12,*) freq_obs(i,j), real(predicted_data(j,i)), imag(predicted_data(j,i))
			end do	
			close(12)
			
		end do

		do i=1,nwave
			sta = stations(i)
			lw=lofw(sta)
			write(filename, '("./rfi_files/FINAL/",a,a13,a1)') 
     &				sta(1:lw),'_observed.asc',char(0)
			open(12, file=filename, status="replace")
			do j=1,nlines(i)
				write(12,*) freq_obs(i,j), real(observed_data(j,i)), imag(observed_data(j,i))
			end do	
			close(12)

			
		end do


	end if





	return
	end








c ----------------------------------------------------------------------------------------------------
c   CALC vc vector
c-----------------------------------------------------------------------------------------------------
	subroutine calc_vc(velocity_abs, beta, strike, dip, vc)
	
	real velocity_abs, beta, strike, dip, theta, pi
	real v0(3), vf(3), vc(3)
	real R(3,3), Rd(3,3), Rs(3,3), I(3,3), S(3,3)

	pi = 4.*atan(1.)
	v0(1) = velocity_abs
	v0(2) = 0.
	v0(3) = 0. 
	R = transpose(reshape((/cos(beta),sin(beta),0.,-sin(beta),cos(beta),0.,0.,0.,1./),shape(R)))
	vf = matmul(v0,R)

	I = transpose(reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),shape(I)))
	Rd = transpose(reshape((/1.,0.,0.,0.,cos(-dip),-sin(-dip),0.,sin(-dip),cos(-dip)/),shape(Rd)))
	theta = pi/2.-strike
	Rs = transpose(reshape((/cos(theta),-sin(theta),0.,sin(theta),cos(theta),0.,0.,0.,1./),shape(Rs)))

	S = matmul(I,Rs)
	S = matmul(S, Rd)
	vc = matmul(vf, transpose(S))

	return
	end

c ----------------------------------------------------------------------------------------------------
c   CALC W matrix
c-----------------------------------------------------------------------------------------------------

	subroutine calc_W(Amax, Amin, phi, dip, strike, W)

	real*4 	M(3,3), S(3,3), W(3,3), I(3,3), R1(3,3), R2(3,3), R3(3,3), U(3,3), Sinv(3,3)
	real*4	S1(3,3)
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

c	write(*,*) R1(1,1), R1(1,2), R1(1,3)
c	write(*,*) R1(2,1), R1(2,2), R1(2,3)
c	write(*,*) R1(3,1), R1(3,2), R1(3,3)
c	write(*,*) "--------------------"

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

c	write(*,*) R2(1,1), R2(1,2), R2(1,3)
c	write(*,*) R2(2,1), R2(2,2), R2(2,3)
c	write(*,*) R2(3,1), R2(3,2), R2(3,3)
c	write(*,*) "--------------------"
c
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

c	write(*,*) R3(1,1), R3(1,2), R3(1,3)
c	write(*,*) R3(2,1), R3(2,2), R3(2,3)
c	write(*,*) R3(3,1), R3(3,2), R3(3,3)
c	write(*,*) "--------------------"


	U = matmul(I,R1)

c	write(*,*) U(1,1), U(1,2), U(1,3)
c	write(*,*) U(2,1), U(2,2), U(2,3)
c	write(*,*) U(3,1), U(3,2), U(3,3)
c	write(*,*) "--------------------"


	U = matmul(U,R2)
c	write(*,*) U(1,1), U(1,2), U(1,3)
c	write(*,*) U(2,1), U(2,2), U(2,3)
c	write(*,*) U(3,1), U(3,2), U(3,3)
c	write(*,*) "--------------------"
	U = matmul(U,R3)
c	write(*,*) U(1,1), U(1,2), U(1,3)
c	write(*,*) U(2,1), U(2,2), U(2,3)
c	write(*,*) U(3,1), U(3,2), U(3,3)
c	write(*,*) "--------------------"

	S = U

c	write(*,*) Sinv(1,1), Sinv(1,2), Sinv(1,3)
c	write(*,*) Sinv(2,1), Sinv(2,2), Sinv(2,3)
c	write(*,*) Sinv(3,1), Sinv(3,2), Sinv(3,3)
c	write(*,*) "--------------------"

	W = matmul(S,M)
	call matinv(S,3,Sinv)
c	write(*,*) W(1,1), W(1,2), W(1,3)
c	write(*,*) W(2,1), W(2,2), W(2,3)
c	write(*,*) W(3,1), W(3,2), W(3,3)
c	write(*,*) "--------------------"
	W = matmul(W,Sinv)

c	write(*,*) W(1,1), W(1,2), W(1,3)
c	write(*,*) W(2,1), W(2,2), W(2,3)
c	write(*,*) W(3,1), W(3,2), W(3,3)
c	write(*,*) "--------------------"


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



c=====================================================================
	subroutine v_given_par(Amax, Amin, beta, phi, duration, v)
	real 	Amax, Amin, beta, phi, duration
	real	AB, r, v
	real	pi, DTOR, RTOD


	r = sqrt(Amin**2 + Amax**2)
	
	AB = abs(r * cos(beta - phi - asin( Amin / r    ) ) )
	v = AB / duration



	return
	end
c======================================================================
	















