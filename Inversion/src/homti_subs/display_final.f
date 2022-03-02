!-------------------------------------------------------------------------
!
!    Subroutine display_final - writes out the final model contained in
!                   array rmodel to the file assigned to
!                                  logical unit lu
!
!    Note: Calls no other routines
!
!       Modified from original NA package
!------------------------------------------------------------------------
!
      subroutine display_final(
     1lu, rmodel, moddim, misfitval,
     1q_alpha, q_beta )

         include 'homti_param.inc'


         real*4        rmodel(maxmoddim)

         real*4        misfitval,
     1                     q_alpha(maxlayer),
     1                     q_beta(maxlayer)


!                        Write out final model

         depth=0.0

         write(lu,*)
         write(lu,*)
         write(lu,*) "==============================================="
         write(lu,*) "==============================================="
         write(lu,801)
         write(lu,*)



         n_param_lines=moddim/3

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



         write(lu,*) "Spatial centroid module: ", rc_module, " angle: ",
     1     rc_angle
         write(lu,*) "Temporal centroid: ",delta_t
         write(lu,*) "Duration: ", duration
         write(lu,*) "Velocity abs: ", vel_abs, " velocity angle: ", 
     1    vel_angle
         write(lu,*) "Amax: ", Amax, " Amin: ", Amin, " phi: ", phi
         write(lu,*) "Strike: ", strike, " dip: ", dip
         write(lu,*) "==============================================="
         write(lu,*) "==============================================="
         write(lu,*)



  801    format('  *** Final model ***' )

         return
      end
