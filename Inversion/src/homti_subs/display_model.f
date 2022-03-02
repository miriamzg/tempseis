!-------------------------------------------------------------------------
!
!    Subroutine display_model - writes out the model contained in
!                   array rmodel and misfit value to
!                   the file assigned to logical
!                   unit lu
!
!    Note: Calls no other routines
!
!------------------------------------------------------------------------
!
      subroutine display_model(
     &lu, imod, rmodel, moddim, misfitval )

         include    'homti_param.inc'


         real*4 rmodel(maxmoddim), pi, DTOR, RDOT, 
     1    vel_angle_rad, phi_rad

         real*4 misfitval

!                        Write out model `rmodel'

         write(lu,801) imod, misfitval

         n_param_lines=moddim/3

         rc_x = rmodel(1)
         rc_y = rmodel(n_param_lines+1)
         rc_z = rmodel(2*n_param_lines+1)
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

         vel_angle_rad = vel_angle * DTOR
         phi_rad = phi * DTOR

!    call v_given_par(Amax, Amin, vel_angle_rad, phi_rad, duration, vel_abs)


!    write(lu,811) "-------------"
         write(lu,*) "Spatial centroid: ", rc_x,rc_y,rc_z
         write(lu,*) "Temporal centroid: ",delta_t
         write(lu,*) "Duration: ", duration
         write(lu,*) "Velocity abs: ", vel_abs, " velocity angle: ",
     1     vel_angle
         write(lu,*) "Amax: ", Amax, " Amin: ", Amin, " phi: ", phi
         write(lu,*) "Strike: ", strike, " dip: ", dip
         write(lu,*) "-------------"



  801    format( '  model:',i5,',   misfit value:',ES20.13 )

c
         return
      end




!=====================================================================
!    subroutine v_given_par(Amax, Amin, beta, phi, duration, v)
!    real     Amax, Amin, beta, phi, duration
!    real    AB, r, v
!    real    pi, DTOR, RTOD
!
!
!    r = sqrt(Amin**2 + Amax**2)
!
!    AB = abs(r * cos(beta - phi - asin( Amin / r    ) ) )
!    v = AB / duration
!
!
!
!    return
!    end
!======================================================================














