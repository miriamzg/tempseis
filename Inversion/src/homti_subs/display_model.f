c-------------------------------------------------------------------------
c
c	Subroutine display_model - writes out the model contained in 
c				   array rmodel and misfit value to 
c				   the file assigned to logical 
c				   unit lu
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c
	subroutine display_model(
     &		lu, imod, rmodel, moddim, misfitval )
c
	include	'homti_param.inc'
c
c
	real*4		rmodel(maxmoddim), pi, DTOR, RDOT, vel_angle_rad, phi_rad
c
	real*4		misfitval
c
c						Write out model `rmodel'
c
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
	
c	call v_given_par(Amax, Amin, vel_angle_rad, phi_rad, duration, vel_abs)

	
c	write(lu,811) "-------------"
	write(lu,*) "Spatial centroid: ", rc_x,rc_y,rc_z
	write(lu,*) "Temporal centroid: ",delta_t 
	write(lu,*) "Duration: ", duration
	write(lu,*) "Velocity abs: ", vel_abs, " velocity angle: ", vel_angle
	write(lu,*) "Amax: ", Amax, " Amin: ", Amin, " phi: ", phi
	write(lu,*) "Strike: ", strike, " dip: ", dip 
	write(lu,*) "-------------"



  801	format( '  model:',i5,',   misfit value:',ES20.13 )

c
	return
	end




c=====================================================================
c	subroutine v_given_par(Amax, Amin, beta, phi, duration, v)
c	real 	Amax, Amin, beta, phi, duration
c	real	AB, r, v
c	real	pi, DTOR, RTOD
c
c
c	r = sqrt(Amin**2 + Amax**2)
c	
c	AB = abs(r * cos(beta - phi - asin( Amin / r    ) ) )
c	v = AB / duration
c
c
c
c	return
c	end
c======================================================================














