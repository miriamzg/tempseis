c-------------------------------------------------------------------------
c
c	Subroutine display_final - writes out the final model contained in 
c				   array rmodel to the file assigned to 
c                                  logical unit lu
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c
	subroutine display_final(
     &		lu, rmodel, moddim, misfitval,  
     &		q_alpha, q_beta )
c
	include 'homti_param.inc'
c
c
	real*4		rmodel(maxmoddim)
c
	real*4		misfitval,
     &                  q_alpha(maxlayer),
     &                  q_beta(maxlayer)
c
c
c						Write out final model
c
	depth=0.0

	write(lu,*)
	write(lu,*)
	write(lu,*) "======================================================="
	write(lu,*) "======================================================="
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
	
	

	write(lu,*) "Spatial centroid module: ", rc_module, " angle: ", rc_angle
	write(lu,*) "Temporal centroid: ",delta_t 
	write(lu,*) "Duration: ", duration
	write(lu,*) "Velocity abs: ", vel_abs, " velocity angle: ", vel_angle
	write(lu,*) "Amax: ", Amax, " Amin: ", Amin, " phi: ", phi
	write(lu,*) "Strike: ", strike, " dip: ", dip 
	write(lu,*) "======================================================="
	write(lu,*) "======================================================="
	write(lu,*)


c
  801	format('  *** Final model ***' )
c
	return
	end
