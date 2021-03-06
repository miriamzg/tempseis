c
c ----------------------------------------------------------------------------
c
c	An example of a driver program that calls the 
c	Neighbourhood algorithm routines to search a 
c	parameter space. See HTML manual for more details.
c					
c                                               M. Sambridge, RSES, ANU.
c						Last updated Sept. 1999.
c
c ----------------------------------------------------------------------------
c
	Program rfi_na

c					MEMORY AND ARRAY SIZES
c
c					The NA routines use the
c					following include file to define
c					all parameters controlling
c					memory required by the arrays.
c
c					A description of each parameter,
c					indicating which ones should be 
c					changed to fit your application
c					can be found in 'na_param.inc' 
c
c	include 	'na_param.inc'
c
c					INPUT AND OUTPUT
c       Input files:
c                   na.in       	Contains options for Neighbourhood 
c					algorithm (see HTML manual)
c
c       Output files:
c                   na.sum      	summary of results
c
c                   na.nad      	NAD (NA direct access file of models)
c
c                   sobol.coeff 	initializing data used for
c                               	quasi-random sequences
c					(This is output for reference only)
c
c       Comments:
c                The NAD file is a direct access compact format file
c                containing all models generated by the neighbourhood
c                algorithm. (NAD files can be read in by multi-dimensional
c                integration program NA-Bayes and plot program S-plot.)
c
c		 Logical units 30-40 are reserved for use by NA subroutines
c		 The user specific routines should not use these values
c		 for logical units.
c		 The NA routines also write to LU 6 as standard out.
c
c ----------------------------------------------------------------------------
c
c					Call NA routine to do the work		
	call na

	stop
	end
c
c ----------------------------------------------------------------------------
c
c	THE FOLLOWING IS AN EXAMPLE OF THE THREE ROUTINES,
c	CALLED BY THE NA ROUTINES, WHICH THE USER MUST SUPPLY. 
c
c	user_init - defines parametrization and can be
c		    used to perform any user initialization 
c		    (e.g. reading in data etc.) 
c
c	forward   - performs forward modelling on an input model
c		    and returns the misfit value of that model.
c
c	writemodels - a routine called after the NA has completed.
c		      All models and misfits are supplied 
c		      to the routine and here the user may write  
c		      then out in their own format.
c		      (Note that the NA routines only write the
c		       models out in a direct access NAD file.
c		       see HTML manual for details.)
c
c	THIS EXAMPLE IS FOR A CRUSTAL RECEIVER FUNCTION INVERSION PROBLEM. 
c
c		In this case the routines perform the following tasks:
c		read in data, define the parametrization,
c		calculate synthetics, define a misfit function and
c		write out various information to ascii files.
c
c	See HTMl manual for more details. An explanation of
c	calling the sequence and meaning is attached to each routine. 
c
c
c       This code is based on work by T. Shibutani (RCEP, DPRI, KYOTO UNIV.)
c	who wrote the receiver function routines.
c
c	Modifications by M. Sambridge (RSES, ANU) for use 
c	with a Neighbourhood algorithm.
c
c						M. Sambridge (RSES, ANU)
c
c-------------------------------------------------------------------------
c
c	user_init - performs all user specific initialization tasks. 
c		    In this case for receiver receiver function inversion.
c
c	Input: - iproc		:processor id (for MPI compilation).
c                                This allows the root node (id=0) to 
c				 emit debug iformation during user_init. 
c				 For serial option id=0.
c
c	Output:
c	      nd		:Number of dimensions in parameter space
c	      ranges(2,nd)	:bounds on parameter space
c	      scales(nd+1)	:scale factors in parameter space
c
c	Comments:
c		 The scale factors are used to transform the model
c		 parameter space into a non-dimensional space. They 
c		 must be provided by the user to put each parameter
c		 on an equal footing.  In the transformation
c		 the i-th parameter becomes,
c		
c		 xi_after = (xi_before-ranges(1,i))/scales(i+1).
c
c		 Note that scales(i+1) contains the scale factor for 
c		 the i-th parameter. The value of scales(1) can be used 
c		 as a shortcut to indicate common choices for scale 
c		 factors. If scales(1) = 0 then all scale factors are
c		 reset to 1; if scales(1) = -1 then the scale factors 
c		 are reset to the parameter range, i.e.,
c
c		 scales(i+1) = ranges(2,i)-ranges(1,i)
c
c		 and for any other value the scale factors are left unchanged.
c		 Note that if scales(1) is set to 0 or -1 then there is
c		 no need to fill up the scales array.
c
c		 This routine is only set up for a 6 layer 24 dimensional 
c		 receiver function problem. 
c
c				M. Sambridge, RSES (ANU), April 1998. 
c
c-------------------------------------------------------------------------
c
	subroutine user_init(nd,ranges,scales)

c						initialize receiver
c						function forward modelling
        include 'rfi_subs/rfi_param.inc'
c
	real*4		ranges(2,*),
     &			scales(*)
c
	real*4		observed_data_r(maxdata,maxwave),
     &                  observed_data_i(maxdata,maxwave),
     &			freq_obs(maxwave, maxdata),
     &                  predicted_data_r(maxdata,maxwave),
     &                  predicted_data_i(maxdata,maxwave),
     &			incident_angle(maxwave),
     &			constant_a(maxwave),
     &			constant_c(maxwave),
     &			time_shift(maxwave),
     &			weight(maxwave),
     &			time_begin(maxwave),
     &			time_end(maxwave)

	complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave)

	real*4		amp_obs_real(maxwave),
     &			amp_obs_imag(maxwave)

c
	integer		ndata(maxwave),
     &			nwave, nd

	logical		verbose,debug,timing,summary
	logical		lroot
c
	character*40	chars,
     &			kname
	character*40    fname(maxwave)

     	character*30	sname(maxwave)
	character*10	stations(maxwave)

	real*4		dS1_freq(3,maxdata)
	real*4		dS1_real(3,maxdata)
	real*4		dS1_imag(3,maxdata)
	complex*8	dS1(maxwave, 3, maxdata)
	real*4		dS2_freq(3,3,maxdata)
	real*4		dS2_real(3,3,maxdata)
	real*4		dS2_imag(3,3,maxdata)
	complex*8	dS2(maxwave, 3,3,maxdata)
	character*100	filename, sta
	real*4		freq_ps(maxwave, maxdata), point_source_real(maxdata), point_source_imag(maxdata)
	complex*8	point_source(maxwave, maxdata)
	integer		nlines(maxwave)


c						Info and Logical unit common 
c						blocks used by NA routines
c
        common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &                lu_nad,verbose,debug,timing,summary

        common /NAMPI/iproc,nproc,lroot

	common /rfi_com/observed_data, weight, 
     &                  time_begin, time_end,  ndata,
     &                  nwave, lu_mod, fname, sname, stations,
     &			dS1, dS2, point_source, freq_obs, nlines,
     &			amp_obs_real, amp_obs_imag



c
c		Set up logical units 
c						LU's for standard input and 
c						output
        lu_rfi  = 15
        lu_out = 6
c
c
c		Set up logical units for rfi files
c
c						LU for input of velocity 
c						model
	lu_vel = 11
c						LU for output of model 
c						parameters
	lu_mod = 65
c


	open(123,file="station_mft.asc", status="unknown")
	open(111,file="best_model.asc", status="unknown")



c						Open rfi files
	open(lu_rfi,file='rfi_files/rfi.in',status='old')
	read(lu_rfi,*)
	read(lu_rfi,*)
	read(lu_rfi,*)

	if(lroot)write(lu_out,*)
	if(lroot)write(lu_out,*)' User routines output'  
	if(lroot)write(lu_out,*)
	if(lroot)write(lu_out,*)' Opening rfi files...'
	read(lu_rfi,'(a)') chars
	lw=lofw(chars)
        write(kname,'("rfi_files/NA_MDL/",a,a1)') chars(1:lw),char(0)
	if(lroot)write(lu_out,*)
	if(lroot)write(lu_out,*) '* Now open ... ',kname(1:lw+17)
	open(lu_vel,file=kname,status='old')
c
	read(lu_rfi,'(a)') chars
	lw=lofw(chars)
        write(kname,'("rfi_files/NA_MDL/",a,a1)') chars(1:lw),char(0)
	if(lroot)write(lu_out,*)
	if(lroot)write(lu_out,*) '* Now open ... ',kname(1:lw+17)
	open(lu_mod, file=kname, status='replace')



c
c						Read in velocity model 
c						and ranges of the parameters
	call read_vmodelNA
     &	     (lu_vel, ranges, scales, nd ) 
	close(lu_vel)



c
c
c						Read in ORFs in SAC format
	read(lu_rfi,*) nwave
	do iw=1,nwave
	  read(lu_rfi,'(a)') chars
	  lw=lofw(chars)
	  stations(iw) = chars(1:lw)
	  write(fname(iw),'("rfi_files/OBS/",a,a3,a1)') chars(1:lw),'_ff',char(0)
	  write(sname(iw),'("rfi_files/SYNT/",a,a3,a1)') chars(1:lw),'_ff',char(0)

	  if(lroot)write(lu_out,*)
	  if(lroot)write(lu_out,*) '* Read in data from ',fname(iw)(1:lw+18)
c
	  read(lu_rfi,*) weight(iw)
	  if(lroot)write(lu_out,*) '* Weight =',weight(iw)
	  if(lroot)write(lu_out,*)' '
	


	  open(666, file=fname(iw), status="old")
c	  write(*,*) fname(iw)
	  read(666,*) npoints
	  read(666,*)
	  do i=1,npoints
		read(666,*) freq_obs(iw,i), observed_data_r(i,iw), observed_data_i(i,iw)
		observed_data(i,iw) = complex(observed_data_r(i,iw),observed_data_i(i,iw))
	  end do
	  close(666)
	  
c	  amp_obs_real(iw) = maxval(abs(real(observed_data(:,iw))))
c	  amp_obs_imag(iw) = maxval(abs(imag(observed_data(:,iw))))


	  amp_obs_real(iw) = 0.0
	  amp_obs_imag(iw) = 0.0
	  do i=1, npoints
		amp_obs_real(iw) = amp_obs_real(iw) + real(observed_data(i,iw))**2
	  	amp_obs_imag(iw) = amp_obs_imag(iw) + imag(observed_data(i,iw))**2

	  end do
c	write(*,*) "station ", stations(iw), amp_obs_real(iw), amp_obs_imag(iw), amp_obs_real(iw) + amp_obs_imag(iw)

	end do


	

		
	close(lu_rfi)


c--------------- read kernels
	do i=1,nwave
		sta = stations(i)
		lw=lofw(sta)
c----------------------------------------------
c		load kernels first order
c----------------------------------------------
		write(filename, '("../src/rfi_subs/kernels/",a,a5,a1)') sta(1:lw),'_dSdx',char(0)
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS1_freq(1,j), dS1_real(1,j), dS1_imag(1,j)
			dS1(i, 1,j)=complex(dS1_real(1,j),dS1_imag(1,j))
		end do
		close(11)

		write(filename, '("../src/rfi_subs/kernels/",a,a5,a1)') sta(1:lw),'_dSdy',char(0)
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS1_freq(2,j), dS1_real(2,j), dS1_imag(2,j)
			dS1(i, 2,j)=complex(dS1_real(2,j),dS1_imag(2,j))
		end do
		close(11)
		write(filename, '("../src/rfi_subs/kernels/",a,a5,a1)') sta(1:lw),'_dSdz',char(0)
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS1_freq(3,j), dS1_real(3,j), dS1_imag(3,j)
			dS1(i, 3,j)=complex(dS1_real(3,j),dS1_imag(3,j))
		end do
		close(11)

c----------------------------------------------
c		load kernels second order
c----------------------------------------------

		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dS2dx2',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(1,1,j), dS2_real(1,1,j), dS2_imag(1,1,j)
			dS2(i, 1,1,j)=complex(dS2_real(1,1,j),dS2_imag(1,1,j))
			
		end do
		close(11)

		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdxdy',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(1,2,j), dS2_real(1,2,j), dS2_imag(1,2,j)
			dS2(i, 1,2,j)=complex(dS2_real(1,2,j),dS2_imag(1,2,j))
		end do
		close(11)

		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdxdz',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(1,3,j), dS2_real(1,3,j), dS2_imag(1,3,j)
			dS2(i, 1,3,j)=complex(dS2_real(1,3,j),dS2_imag(1,3,j))
		end do
		close(11)

		dS2_freq(2,1,:) = dS2_freq(1,2,:)
		dS2(i, 2,1,:) = dS2(i, 1,2,:)
	
		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dS2dy2',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(2,2,j), dS2_real(2,2,j), dS2_imag(2,2,j)
			dS2(i, 2,2,j)=complex(dS2_real(2,2,j),dS2_imag(2,2,j))
		end do
		close(11)

		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdydz',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(2,3,j), dS2_real(2,3,j), dS2_imag(2,3,j)
			dS2(i, 2,3,j)=complex(dS2_real(2,3,j),dS2_imag(2,3,j))
			
		end do
		close(11)

		dS2_freq(3,1,:) = dS2_freq(1,3,:)
		dS2(i, 3,1,:) = dS2(i, 1,3,:)

		dS2_freq(3,2,:) = dS2_freq(2,3,:)
		dS2(i, 3,2,:) = dS2(i, 2,3,:)


		write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dS2dz2',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11, *) nlines(i)
		read(11,*)
		do j=1,nlines(i)
			read(11,*) dS2_freq(3,3,j), dS2_real(3,3,j), dS2_imag(3,3,j)
			dS2(i, 3,3,j)=complex(dS2_real(3,3,j),dS2_imag(3,3,j))
		end do
		close(11)
c-----------------------------------------------------------------------------------
c		read point source solution
c-----------------------------------------------------------------------------------

		write(filename, '("../src/rfi_subs/point_source/",a,a3,a1)') sta(1:lw),'_ps',char(0)
c		write(*,*) filename
		open(11, file=filename, status="old")
		read(11,*) nlines(i)
		read(11,*)
		do j=1, nlines(i)
			read(11,*) freq_ps(i, j), point_source_real(j), point_source_imag(j)
			point_source(i, j)=complex(point_source_real(j), point_source_imag(j))
		end do
		close(11)

	end do







c-------------------------------
c
        if(lroot)write(lu_out,*)
c	if(lroot)write(lu_out,101)

 100    format(1x,60("-")//)
 101    format(/1x,60("-")/)




	return
	end
c
c-------------------------------------------------------------------------
c
c	forward - performs forward modelling for user supplied problem.
c		  In this case it calculates predicted receiver function
c		  for a single model and calculates the misfit measure
c		  between observation and prediction.
c
c	Input: 
c	      nd		:Number of dimensions in parameter space
c	      model(nd)		:input velocity model
c
c	Output:
c	      lppd		:negative log(ppd)
c
c	Comments:
c		 This routine allows the user to perform the forward 
c		 modelling and define an a posterior probability density
c		 function using the resulting mismatch with the observed
c		 data. Note that the output quanity is the negative
c		 natural log of the ppd. This is equivalent to the 
c		 misfit function in most problems. 
c
c				M. Sambridge, RSES (ANU), April 1998. 
c
c-------------------------------------------------------------------------
c
	subroutine forward(nd,model,lppd)

c						initialize receiver
c						function forward modelling
        include 'rfi_subs/rfi_param.inc'

	real*4		lppd, misfitval

	real*4		model(nd)

	real*8		rmodel(max_nd)
	complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave)
c
	real*4		amp_obs_real(maxwave),
     &			amp_obs_imag(maxwave)

	real*4		observed_data_r(maxdata,maxwave),
     &                  observed_data_i(maxdata,maxwave),
     &                  predicted_data_r(maxdata,maxwave),
     &                  predicted_data_i(maxdata,maxwave),
     &			freq_obs(maxwave, maxdata),
     &			weight(maxwave),
     &			time_begin(maxwave),
     &			time_end(maxwave),
     &			Scalc_real(maxdata), Scalc_imag(maxdata)
c
	integer		ndata(maxwave), s
	integer		nwave, nd
c	character*50	filename

	logical		verbose,debug,timing,summary

	character*40	fname(maxwave)
     	character*30	sname(maxwave)
	character*300   command
	character*10	stations(maxwave)

	real*4		dS1_freq(3,maxdata)
	real*4		dS1_real(3,maxdata)
	real*4		dS1_imag(3,maxdata)
	complex*8	dS1(maxwave, 3, maxdata)
	real*4		dS2_freq(3,3,maxdata)
	real*4		dS2_real(3,3,maxdata)
	real*4		dS2_imag(3,3,maxdata)
	complex*8	dS2(maxwave, 3,3,maxdata)
	character*100	filename, sta
	real*4		freq_ps(maxwave, maxdata), point_source_real(maxdata), point_source_imag(maxdata)
	complex*8	point_source(maxwave, maxdata)
	integer		nlines(maxwave)

c
c						Info and Logical unit common 
c						blocks used by NA routines
c
        common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &                lu_nad,verbose,debug,timing,summary


	common /rfi_com/observed_data, weight, 
     &                  time_begin, time_end,  ndata,
     &                  nwave, lu_mod, fname, sname, stations,
     &			dS1, dS2, point_source, freq_obs, nlines,
     &			amp_obs_real, amp_obs_imag




c
c		Set up logical units 
c						LU's for standard input and 
c						output
        do j=1,nd
           rmodel(j) = dble(model(j))
        end do
c                                               Perform forward modelling
c                                               on model `rmodel'.



c	write(*,*) "NDATA", ndata(1)


        call forward_modelling(
     &      	rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &		dS1, dS2, point_source, freq_obs, nlines, .false., .false.)

	
c	write(*,*) "NDATA", ndata(1)

c	Debug: write down synthetic file and plot it
c	do s=1,nwave
c		open(91, file=sname(s), status="replace")
c		write(91,*) ndata(s)
c		write(91,*) "freq		real			imaginary"
c		do i=1,ndata(s)
c			write(91,*) freq(i,s), real(predicted_data(i,s)), imag(predicted_data(i,s))
c		end do
c		close(91)
c	end do
c	do s=1,nwave
c		write(filename,'(a20,a4,a1)') trim(sname(s)),"_obs",char(0)
c		write(*,*) "debug ", filename
c		open(71, file=filename, status="replace")
c		write(71,*) ndata(s)
c		write(71,*) "freq		real			imaginary"
c		do i=1,ndata(s)
c			write(71,*) freq(i,s), real(observed_data(i,s)), imag(observed_data(i,s))
c		end do
c		close(71)
c
c		command = "python plot_data_synt.py " // sname(s) // " " // filename
c		call system(command)
c		
c	end do
c	stop

c                                               Calculate misfit function


        call calcmisfit(
     &       predicted_data, 
     &       observed_data, nlines,
     &       weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .false.)
c
	
	lppd =  misfitval * 1



	
	write(*,*) "Misfit: ", lppd
	write(123,*) "ooooooooooooooooooooo"
	write(123,*) "Tot Misfit: ", lppd
c	stop
c	write(*,*) "----------------------"
	write(123,*) "======================================"
	

	return
	end
c
c-------------------------------------------------------------------------
c
c	writemodels - user supplied routine to write out models produced
c		      by Neighbourhood algorithm in user's own format.
c
c	Input: 
c	      nd		  :number of dimensions in parameter space
c	      ntot		  :number of models generated by NA 
c	      models(nd,ntot)     :models generated by NA
c	      misfit 		  :array of model misfits (-lppd's)
c	      ns1  		  :initial sample size used by NA
c	      ns2  		  :normal sample size used by NA
c	      itmax  		  :number of iterations
c	      nh_max  		  :maximum length of nad file header 
c
c	Output:         
c	      nh  		  :length of nad file header 
c	      header(nh)  	  :character string containing nad file header
c
c	Comments:
c	         This routine allows the user to write out the models
c		 produced by the NA in their own format. The three 
c		 parameters ns1, ns2 and itmax define the samples sizes
c		 and number of iterations taken by the NA. Initially ns1
c		 random samples are generated, and then itmax sets of ns2
c		 models are generated and so ntot = ns1 + ns2*itmax. 
c
c		 The user may if they wish write their own information 
c		 into the `user reserved' part of NAD file header.
c		 (See example below for details).
c
c				M. Sambridge, RSES (ANU), Oct 1999. 
c
c-------------------------------------------------------------------------
c
	subroutine writemodels
     &             (nd, ntot, models, misfit, ns1, ns2, itmax,
     &              nh_max, nh, header)

        include 'rfi_subs/rfi_param.inc'

c						NA variables and arrays
	real*4		models(nd,*)
	real*4		misfit(ntot)
	real*4		mfitmin
	real*4		mfitminc
	real*4		mfitmean
	real*4		lppd, misfitval

        character*(*)   header

	logical		verbose,debug,timing,summary

c						Receiver function variables
c						and arrays
	real*8		rmodel(max_nd)

	complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave)

	real*4		amp_obs_real(maxwave),
     &			amp_obs_imag(maxwave)

	real*4		weight(maxwave),
     &			time_begin(maxwave),
     &			time_end(maxwave),
     &			freq_obs(maxwave, maxdata)

	integer		ndata(maxwave), nd

	character*40    fname(maxwave)
     	character*30	sname(maxwave)
	character*10	stations(maxwave)

	
	complex*8	dS1(maxwave, 3, maxdata)
	complex*8	dS2(maxwave, 3,3,maxdata)
	character*100	filename, sta
	real*4		freq_ps(maxdata), point_source_real(maxdata), point_source_imag(maxdata)
	complex*8	point_source(maxwave, maxdata)
	integer		nlines(maxwave)

c
c						Info and Logical unit common 
c						block used by NA routines
c
        common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &                lu_nad,verbose,debug,timing,summary


	common /rfi_com/observed_data, weight, 
     &                  time_begin, time_end,  ndata,
     &                  nwave, lu_mod, fname, sname, stations,
     &			dS1, dS2, point_source, freq_obs, nlines,
     &			amp_obs_real, amp_obs_imag




	nlayers = 6
c						write out models
c						at each iteration
	mfitmin = misfit(1)
	ns = ns1
        np = 0
	mopt = 1
c						turn off writing to 
c						standard out by setting
c						lu to zero
	lu_out2 = lu_out
	lu_out2 = 0
	write(lu_mod,*)ns1,' Number of samples in starting pool'
	write(lu_mod,*)ns2,' Number of new samples per iteration'
	write(lu_mod,*)itmax,' Number of iterations'
	write(*,*) "lu_mod: ", lu_mod
c
c						loop over iterations
	do it=1,itmax+1
	   mfitminc = misfit(np+1)
	   mfitmean = 0.0
c						find minimum and 
c						mean misfit
           do i=1,ns
              jj = np + i
	      if(misfit(jj).lt.mfitmin)then
                 mfitmin = misfit(jj)
                 mopt = jj
              end if
              mfitminc = min(mfitminc,misfit(jj))
              mfitmean = mfitmean + misfit(jj)
           end do
	   mfitmean = mfitmean/ns
c						write out all models
c						from this iteration
c						to receiver function
c						output file.

           write(lu_mod,801) it-1, mfitmin, mfitmean, mfitminc
           do i=1,ns
                jj = np + i
                call display_model
     &               (lu_mod, i, models(1,jj), nd, misfit(jj) )
           end do
	   np = np + ns
           ns = ns2
c                                               Write out summary of
c                                               optimum model to
c						receiver function file
c						and NA summary file.
c
           call output_summary(
     &                  lu_out2, lu_sum, it-1, models(1,mopt), nd,
     &                  np, mfitmin, mfitmean, mfitminc, mopt)
        end do
c                                               Write out final model
c
        call display_final(
     &          lu_out, models(1,mopt), nd, mfitmin )
c
        call display_final(
     &          lu_sum, models(1,mopt), nd, mfitmin )
c
        call display_final(
     &          lu_mod, models(1,mopt), nd, mfitmin)
c
        do j=1,nd
           rmodel(j) = dble(models(j,mopt))
        end do



	write(*,*) "PEPPAAAAAAAAAAA"
	write(111,*) "************************  INITIAL MISFIT ************************ "
	call forward_modelling(
     &      	rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &		dS1, dS2, point_source, freq_obs, nlines, .false., .true.)


	call calcmisfit(
     &       predicted_data, 
     &       observed_data, nlines,
     &       weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .true.)
	
	write(111,*) "***"
	write(111,*) "Total misfit: ", misfitval
	write(111,*) "***"
	


c						repeat forward modelling 
c						for optimum model


	
	write(111,*) "************************ BEST MODEL ************************ "
        call forward_modelling(
     &      	rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &		dS1, dS2, point_source, freq_obs, nlines, .true., .false.)


        call calcmisfit(
     &       predicted_data, 
     &       observed_data, nlines,
     &       weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .true.)
c
	
	

	write(111,*) "***"
	write(111,*) "total misfit: ", misfitval
	write(111,*) "***"
	close(111)





c						write out predicted 
c						receiver function 
c						for optimum model
c        do iw=1,nwave
c          call writedata(
c     &          fname(iw), ndata(iw),
c     &          observed_data(1,iw), predicted_data(1,iw),
c     &          time_begin(iw),
c     &          constant_a(iw), constant_c(iw), 1./fs )
c        end do

c
c						write RFI component
c						of header for NAD file
c						(can leave blank
c						 and set nh to 0)
c	nh = 0
c
c						calculate size of RFI header
        do i=1,nlayers
           k1 = 2*(i-1)*6 + 25
           k2 = k1 + 11
        end do
        nh = k2
 
        if(nh.gt.nh_max)then
           write(*,*)
           write(*,*)' Error - header array too small'
           write(*,*)
           write(*,*)'         current size = ',nh_max
           write(*,*)'        required size = ',nh
           write(*,*)
           stop
        end if
c						write info into 
c						character string

	write(header(1:24),fmt='(4i6)')ns1,ns2,itmax,nlayers
 
        do i=1,nlayers
           k1 = 2*(i-1)*6 + 25
           k2 = k1 + 11
c           write(header(k1:k2),fmt='(2(1x,f5.0))')
c    &           q_alpha(i),q_beta(i)
        end do
 
 
  801      format( 'iteration:',i5,',  misfit: min=',ES10.3,
     &             ', mean=',ES10.3,', minc=',ES10.3 )

	return
	end
