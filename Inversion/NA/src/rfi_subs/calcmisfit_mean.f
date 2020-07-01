c-------------------------------------------------------------------------
c
c	Subroutine calcmisfit - calculates a misfit value between 
c				observed data and predicted data
c				plus model roughness
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c
	subroutine calcmisfit( predicted_data, observed_data, ndata,
     &       weight, nwaves, totmisfit, fname, amp_obs_real, amp_obs_imag )
c
c
	include 'rfi_param.inc'
c
c
	complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave),
     &			observed_phase,
     &			predicted_phase

	real*4		weight(maxwave),
     &                  abs_value_obs, abs_value_synt 
	real*4		real_mft, imag_mft
     	real*4		totmisfit, misfitval, aval_real, aval_imag, aval, mft_abs, mft_fft
	real*4		mft_phase, mft_absolute, amplitude
	real*4		totmisfit_P,totmisfit_S,totmisfit_W
	integer		n_P, n_S, n_W
c
	integer         ndata(maxwave), nwaves
	character*40	fname(maxwave)
	character*40	instring, string_r, string1, string2, delim
	character	wavetype

	real*4		amp_obs_real(maxwave),
     &			amp_obs_imag(maxwave)

	pi =4.D0*DATAN(1.D0)
c
c						misfit between observed and
c						predicted
c	misfitval=0.0
	totmisfit=0.0

	totmisfit_P = 0.0
	totmisfit_S = 0.0
	totmisfit_W = 0.0

	n_P = 0
	n_S = 0
	n_W = 0

	do iw=1,nwaves
		instring = fname(iw)
		call reverse_string(instring, string_r)
		delim = "_"
		call split_string(string_r, string1, string2, delim)
		wavetype = string2(1:1)
		


	

		aval=0.0
		aval_real = 0.0
		aval_imag = 0.0
		mft_fft = 0.0
		mft_abs = 0.0
		mft_phase = 0.0
		mft_absolute = 0.0
c		amplitude = (max_obs_real(iw) + max_obs_imag(iw))/2.


c		amplitude = 0.0
c		do i=1,ndata(iw)
c			amplitude = amplitude +  real(observed_data(i,iw))**2 + imag(observed_data(i,iw))**2
c		end do

c		write(*,*) "amp", amplitude, 
		

		do i=1,ndata(iw)
			
c			----------------------------
c			Original code, do not edit
			real_mft = (real(predicted_data(i,iw)) - real(observed_data(i,iw)))**2 
			imag_mft = (imag(predicted_data(i,iw)) - imag(observed_data(i,iw)))**2 
c			aval = aval + (real_mft + imag_mft) * weight(iw)
c			----------------------------
c			aval = aval + (real_mft/max_obs_real(iw) + imag_mft/max_obs_imag(iw))
			aval = aval + (real_mft + imag_mft) 
c			aval_real = aval_real + real_mft
c			aval_imag = aval_imag + imag_mft



		end do


c		aval_real = aval_real / amp_obs_real(iw)
c		aval_imag = aval_imag / amp_obs_imag(iw)

c		aval = aval_real + aval_imag


		aval = aval / (amp_obs_real(iw) + amp_obs_imag(iw))
		
c		write(*,*) "STATION: ", fname(iw)(15:25) , " Misfit: ", aval
		write(123,*) "STATION: ", fname(iw)(15:26) , " Misfit: ", aval


c		totmisfit = totmisfit + (0.97 * mft_abs + 0.3 * mft_fft)
c		totmisfit = totmisfit + mft_abs		

c		*** ORGINAL LINE ***
c		totmisfit = totmisfit + aval

		if (wavetype .eq. "P") then
			totmisfit_P = totmisfit_P + aval
			n_P = n_P + 1
		end if
		if (wavetype .eq. "S") then
			totmisfit_S = totmisfit_S + aval
			n_S = n_S + 1
		end if
		if (wavetype .eq. "W") then
			totmisfit_W = totmisfit_W + aval
			n_W = n_W + 1
		end if

		if (n_P .eq. 0) then
			n_P = 1
		end if	
		if (n_S .eq. 0) then
			n_S = 1
		end if	
		if (n_W .eq. 0) then
			n_W = 1
		end if		

		totmisfit = (totmisfit_P/n_P) + (totmisfit_S/n_S) + (totmisfit_W/n_W) 

c		*****************************

c		totmisfit = totmisfit + mft_absolute + mft_phase
		
		
	end do


	write(*,*) "mft P: ", totmisfit_P, " mft S: ", totmisfit_S, " mft W: ", totmisfit_W
	write(*,*) "n P: ", n_P, " n S: ", n_S, " n W: ", n_W
	write(*,*) "mft P scaled : ", totmisfit_P/n_P, " mft S scaled : ", totmisfit_S/n_S, " mft W scaled : ", totmisfit_W/n_W

c	close(666)

	return
	end


	subroutine reverse_string(instring, outstring)
	character*40	instring, outstring
	character	temp
	integer		i, length

	
	length = len_trim(instring) ! ignores trailing blanks. 
                            ! use len(string) to reverse those as well

	outstring = instring
  	do i = 1, length/2
     		temp = outstring(i:i)
     		outstring(i:i) = outstring(length+1-i:length+1-i)
    		outstring(length+1-i:length+1-i) = temp
  	end do

	
	end subroutine reverse_string



  ! split a string into 2 either side of a delimiter token
	SUBROUTINE split_string(instring, string1, string2, delim)
	    CHARACTER(40) :: instring,delim
	    CHARACTER(40),INTENT(OUT):: string1,string2
	    INTEGER :: index

	    instring = adjustl(instring)

	    index = SCAN(instring, delim)
	    string1 = instring(1:index-1)
	    string2 = instring(index+1:)

	END SUBROUTINE split_string













