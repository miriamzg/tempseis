      subroutine forward_modelling(
     1           rmodel, moddim, ndata, nwave, stations,
     1           predicted_data, observed_data,
     1           dS1, dS2, point_source, freq_obs, nlines, output_final,
     1           initial_misfit)

         include 'homti_param.inc'

         parameter     ( v60 = 8.043 )
         parameter     ( rad = 0.017453292 )

         integer moddim
         real*8 rmodel(maxmoddim)
         real*4 predicted_data_r(maxdata,maxwave),
     1          predicted_data_i(maxdata,maxwave)
         complex*8 predicted_data(maxdata,maxwave),
     1             observed_data(maxdata,maxwave)
         integer ndata(maxwave), nstations
         real Amax, Amin,  phi, Lmax, Lmin, pi,
     1        DTOR, RTOD, duration, Amin_Amax_ratio
         character*200 command, filename
         character*10 stations(maxwave)
         character*8 sta
         character*50 fault_file
         real rc_lon, rc_lat, rc_depth, dip, strike, rc(3), vc(3)
         real rc_module, rc_angle, rc_angle_rad, vel_abs
         real dip_rad, strike_rad, phi_rad, vel_angle, vel_angle_rad
         character*10 junk

         real a(3,3), y(3,3), W(3,3)
         complex*8 W_imag(3,3)
         real*4 freq_obs(maxwave, maxdata)
         real*4 Scalc_real(maxdata), Scalc_imag(maxdata)
         complex*8 i_unit
         complex*8 AA(maxdata),BB(maxdata),CC(maxdata),
     1             DD(maxdata),EE(maxdata)
         complex*8 tensor_dot, dot, DD_dot
         real*4 tensor_dot_real, tensor_dot_imag
         logical output_final, output_temp, initial_misfit

         complex*8 dS1(maxwave, 3, maxdata)
         complex*8 dS2(maxwave, 3,3,maxdata)
         complex*8 point_source(maxwave, maxdata)
         integer nlines(maxwave)

         i_unit = (0,1)
         output_temp = .false.
         n_param_lines = 6
         rc_module = rmodel(1)
         rc_angle = rmodel(n_param_lines+1)
         delta_t = rmodel(2)
         duration = rmodel(3)
         vel_abs = rmodel(4)
         vel_angle = rmodel(n_param_lines+4)
         Amax = rmodel(5)
         Amin_Amax_ratio = rmodel(n_param_lines+5)
         Amin = Amax * Amin_Amax_ratio

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

         write(*,*) "-------------"
         write(*,*) "Spatial centroid module: ", rc_module,
     1              " angle: ", rc_angle
         write(*,*) "Temporal centroid: ",delta_t
         write(*,*) "Duration: ", duration * sqrt(3.)
         write(*,*) "Velocity abs: ", vel_abs,
     1              " velocity angle: ", vel_angle
         write(*,*) "Amax: ", Amax * sqrt(3.), " Amin: ", 
     1              Amin * sqrt(3.),
     1              " phi: ", phi
         write(*,*) "Strike: ", strike, " dip: ", dip
         write(*,*) "-------------"


         write(123,*) "Spatial centroid module: ", rc_module,
     1                " angle: ", rc_angle
         write(123,*) "Temporal centroid: ",delta_t
         write(123,*) "Duration: ", duration * sqrt(3.)
         write(123,*) "Velocity abs: ", vel_abs,
     1                " velocity angle: ", vel_angle
         write(123,*) "Amax: ", Amax * sqrt(3.),
     1                " Amin: ", Amin * sqrt(3.),
     1                " phi: ", phi
         write(123,*) "Strike: ", strike, " dip: ", dip
         write(123,*) "***********************"


         if (output_final) then
            write(111,*) "Spatial centroid module: ", rc_module,
     1                   " angle: ", rc_angle
            write(111,*) "Temporal centroid: ",delta_t
            write(111,*) "Duration: ", duration * sqrt(3.)
            write(111,*) "Velocity abs: ", vel_abs,
     1                   " velocity angle: ", vel_angle
            write(111,*) "Amax: ", Amax * sqrt(3.),
     1                   " Amin: ", Amin * sqrt(3.),
     1                   " phi: ", phi
            write(111,*) "Strike: ", strike, " dip: ", dip
            write(111,*) "***********************"
         end if


         call calc_vc(vel_abs, vel_angle_rad, strike_rad, dip_rad, vc)

         call calc_vc(rc_module, rc_angle_rad, strike_rad, dip_rad, rc)

         do i=1,nwave
            sta = stations(i)
            lw=lofw(sta)

            call calc_W(Amax, Amin, phi_rad, dip_rad, strike_rad, W)

            do j=1,nlines(i)
               AA(j)= i_unit * freq_obs(i,j) * (-delta_t) *
     1               point_source(i,j)
!------------------------------------------------------------------------
               BB(j) = 0.0
               do k=1,3
                  BB(j) = BB(j) + rc(k) * dS1(i,k,j)
               end do
!------------------------------------------------------------------------
               CC(j) = -(freq_obs(i, j)**2 / 2.0) * (duration/2.0)**2 *
     1                point_source(i,j)
!------------------------------------------------------------------------
               DD_dot = 0.0
               do k=1,3
                  DD_dot = DD_dot +  (vc(k)*dS1(i,k,j))
               end do
               DD(j) = -i_unit * freq_obs(i,j) * (duration/2.0)**2 *
     1                 DD_dot
!------------------------------------------------------------------------
               EE(j) = (0.0,0.0)
               tensor_dot = (0.0,0.0)
               do k=1,3
                  do l=1,3
                     tensor_dot = tensor_dot + W(k,l) * dS2(i,k,l,j)
                  end do
               end do
               EE(j) =  0.5 * tensor_dot
!------------------------------------------------------------------------
               if (initial_misfit) then
                  predicted_data(j,i) = point_source(i,j)
               else
                  predicted_data(j,i) = point_source(i,j) + AA(j) + 
     1                                  BB(j) + CC(j) + DD(j) + EE(j)
               end if
            end do
         end do

         if (output_temp) then
            do i=1,nwave
               sta = stations(i)
               lw=lofw(sta)
               write(filename, '("TEMP/",a,a14,a1)')
     1                        sta(1:lw),'_predicted.asc',char(0)
               open(12, file=filename, status="replace")
               do j=1,nlines(i)
                  write(12,*) freq_obs(i,j),
     1                        real(predicted_data(j,i)),
     1                        imag(predicted_data(j,i))
               end do
               close(12)

            end do

            do i=1,nwave
               sta = stations(i)
               lw=lofw(sta)
               write(filename, '("TEMP/",a,a13,a1)')
     1                        sta(1:lw),'_observed.asc',char(0)
               open(12, file=filename, status="replace")
               do j=1,nlines(i)
                  write(12,*) freq_obs(i,j),
     1                        real(observed_data(j,i)),
     1                        imag(observed_data(j,i))
               end do
               close(12)

            end do
         end if


         if (output_final) then
            do i=1,nwave
               sta = stations(i)
               lw=lofw(sta)
               write(filename, '("FINAL/",a,a19,a1)')
     &         sta(1:lw),'_best_predicted.asc',char(0)
               open(12, file=filename, status="replace")
               do j=1,nlines(i)
                  write(12,*) freq_obs(i,j),
     1                        real(predicted_data(j,i)),
     1                        imag(predicted_data(j,i))
               end do
               close(12)

            end do

            do i=1,nwave
               sta = stations(i)
               lw=lofw(sta)
               write(filename, '("FINAL/",a,a13,a1)')
     &         sta(1:lw),'_observed.asc',char(0)
               open(12, file=filename, status="replace")
               do j=1,nlines(i)
                  write(12,*) freq_obs(i,j),
     1                        real(observed_data(j,i)),
     1                        imag(observed_data(j,i))
               end do
               close(12)
            end do
         end if

         return
      end




!-----------------------------------------------------------------------
!   CALC vc vector
!------------------------------------------------------------------------
      subroutine calc_vc(velocity_abs, beta, strike, dip, vc)

         real velocity_abs, beta, strike, dip, theta, pi
         real v0(3), vf(3), vc(3)
         real R(3,3), Rd(3,3), Rs(3,3), I(3,3), S(3,3)

         pi = 4.*atan(1.)
         v0(1) = velocity_abs
         v0(2) = 0.
         v0(3) = 0.
         R = transpose(reshape((/cos(beta),sin(beta),0.,-sin(beta),
     1                 cos(beta),0.,0.,0.,1./),shape(R)))
         vf = matmul(v0,R)

         I = transpose(reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),shape(I)))
         Rd = transpose(reshape((/1.,0.,0.,0.,cos(-dip),-sin(-dip),0.,
     1                  sin(-dip),cos(-dip)/),shape(Rd)))
         theta = pi/2.-strike
         Rs = transpose(reshape((/cos(theta),-sin(theta),0.,sin(theta),
     1                  cos(theta),0.,0.,0.,1./),shape(Rs)))

         S = matmul(I,Rs)
         S = matmul(S, Rd)
         vc = matmul(vf, transpose(S))

         return
      end

!------------------------------------------------------------------------
!   CALC W matrix
!------------------------------------------------------------------------

      subroutine calc_W(Amax, Amin, phi, dip, strike, W)

         real*4 M(3,3), S(3,3), W(3,3),
     1          I(3,3), R1(3,3), R2(3,3),
     1          R3(3,3), U(3,3), Sinv(3,3)
         real*4 S1(3,3)
         real phi, dip, strike, a, Amax, Amin, Lmax, Lmin, teta

         Lmax = (Amax/2.)**2
         Lmin = (Amin/2.)**2
         pi = 4.*atan(1.)

!   Matrix containing the eigenvalues
         M(1,1) = Lmax
         M(1,2) = 0.0
         M(1,3) = 0.0
         M(2,1) = 0.0
         M(2,2) = Lmin
         M(2,3) = 0.0
         M(3,1) = 0.0
         M(3,2) = 0.0
         M(3,3) = 0.0

!   Identity matrix
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

!   Rotation matrix strike
         R1(1,1) = cos(teta)
         R1(1,2) = -sin(teta)
         R1(1,3) = 0.0
         R1(2,1) = sin(teta)
         R1(2,2) = cos(teta)
         R1(2,3) = 0.0
         R1(3,1) = 0.0
         R1(3,2) = 0.0
         R1(3,3) = 1.0

!   Rotation matrix dip
         R2(1,1) = 1.0
         R2(1,2) = 0.0
         R2(1,3) = 0.0
         R2(2,1) = 0.0
         R2(2,2) = cos(-dip)
         R2(2,3) = -sin(-dip)
         R2(3,1) = 0.0
         R2(3,2) = sin(-dip)
         R2(3,3) = cos(-dip)

!   Rotation matrix phi
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

         W = matmul(S,M)
         call matinv(S,3,Sinv)
         W = matmul(W,Sinv)

         return
      end



!------------------------------------------------------------------------
!
!     purpose: Inverse computation of a matrix
!
!     Input and output varibles:
!             input         n  size of matrix
!             input         a  input matrix
!             output        y  output inverse matrix
!
!     Use:
!             ludcmp        LU docomposition
!             lubksb        find inverse of matrix by columns
!
!     Warning: After matinv(a,n,y) is called, the input matrix
!              'a' would be changed (used as working space)
!              unless it is a diagnal matix. So if 'a' is
!              not a diagnal matrix and will be used after
!              calling matinv.f, 'a' should be saved before
!              calling matinv.f.
!
!------------------------------------------------------------------------

      subroutine matinv(a,n,y)

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

         return
      end

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
 11         continue
            if (aamax.eq.0.) pause 'singular matrix in ludcmp'
            vv(i)=1./aamax
 12      continue
         do 19 j=1,n
            do 14 i=1,j-1
               sum=a(i,j)
               do 13 k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
 13            continue
               a(i,j)=sum
 14         continue
            aamax=0.
            do 16 i=j,n
               sum=a(i,j)
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
 15            continue
               a(i,j)=sum
               dum=vv(i)*abs(sum)
               if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
               endif
 16         continue
            if (j.ne.imax)then
               do 17 k=1,n
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
 17            continue
               d=-d
               vv(imax)=vv(j)
            endif
            indx(j)=imax
            if(a(j,j).eq.0.)a(j,j)=TINY
            if(j.ne.n)then
               dum=1./a(j,j)
               do 18 i=j+1,n
                  a(i,j)=a(i,j)*dum
 18            continue
            endif
 19      continue
         return
      END

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
 11            continue
            else if (sum.ne.0.) then
               ii=i
            endif
            b(i)=sum
 12      continue
         do 14 i=n,1,-1
            sum=b(i)
            do 13 j=i+1,n
               sum=sum-a(i,j)*b(j)
 13         continue
            b(i)=sum/a(i,i)
 14      continue
         return
      END

!-----------------------------------------------------------------------
      subroutine v_given_par(Amax, Amin, beta, phi, duration, v)
         real  Amax, Amin, beta, phi, duration
         real AB, r, v
         real pi, DTOR, RTOD

         r = sqrt(Amin**2 + Amax**2)
         AB = abs(r * cos(beta - phi - asin( Amin / r    ) ) )
         v = AB / duration

         return
      end
!-----------------------------------------------------------------------
















