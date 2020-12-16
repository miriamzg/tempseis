      Program rfi_na

!					Call NA routine to do the work
      call na

      stop
      end
!------------------------------------------------------------------------
      subroutine user_init(nd,ranges,scales)

      include 'rfi_subs/rfi_param.inc'

      real*4 ranges(2,*), scales(*)
      real*4 observed_data_r(maxdata,maxwave),
     1       observed_data_i(maxdata,maxwave),
     1       freq_obs(maxwave, maxdata),
     1       predicted_data_r(maxdata,maxwave),
     1       predicted_data_i(maxdata,maxwave),
     1       incident_angle(maxwave),
     1       constant_a(maxwave),
     1       constant_c(maxwave),
     1       time_shift(maxwave),
     1       weight(maxwave),
     1       time_begin(maxwave),
     1       time_end(maxwave)

      complex*8 observed_data(maxdata,maxwave),
     1           predicted_data(maxdata,maxwave)

      real*4 amp_obs_real(maxwave),amp_obs_imag(maxwave)


      integer ndata(maxwave),nwave,nd

      logical verbose,debug,timing,summary,lroot

      character*40 chars,kname
      character*40 fname(maxwave)
      character*30 sname(maxwave)
      character*10 stations(maxwave)

      real*4 dS1_freq(3,maxdata)
      real*4 dS1_real(3,maxdata)
      real*4 dS1_imag(3,maxdata)
      complex*8 dS1(maxwave, 3, maxdata)
      real*4 dS2_freq(3,3,maxdata)
      real*4 dS2_real(3,3,maxdata)
      real*4 dS2_imag(3,3,maxdata)
      complex*8 dS2(maxwave, 3,3,maxdata)
      character*100 filename, sta
      real*4 freq_ps(maxwave, maxdata), point_source_real(maxdata)
      real*4 point_source_imag(maxdata)
      complex*8 point_source(maxwave, maxdata)
      integer nlines(maxwave)


!						Info and Logical unit common
!						blocks used by NA routines

      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     1              lu_nad,verbose,debug,timing,summary

      common /NAMPI/iproc,nproc,lroot

      common /rfi_com/observed_data, weight,
     1                time_begin, time_end,  ndata,
     1                nwave, lu_mod, fname, sname, stations,
     1                dS1, dS2, point_source, freq_obs, nlines,
     1                amp_obs_real, amp_obs_imag

!
!		Set up logical units
!       LU's for standard input and output
      lu_rfi  = 15
      lu_out = 6

!		LU for input of velocity model
      lu_vel = 11
!		LU for output of model parameters
      lu_mod = 65

      open(123,file="station_mft.asc", status="unknown")
      open(111,file="best_model.asc", status="unknown")

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

      read(lu_rfi,'(a)') chars
      lw=lofw(chars)
      write(kname,'("rfi_files/NA_MDL/",a,a1)') chars(1:lw),char(0)
      if(lroot)write(lu_out,*)
      if(lroot)write(lu_out,*) '* Now open ... ',kname(1:lw+17)
      open(lu_mod, file=kname, status='replace')


!		Read in velocity model and ranges of the parameters
      call read_vmodelNA(lu_vel, ranges, scales, nd)
      close(lu_vel)


!		Read in observed data in SAC format
      read(lu_rfi,*) nwave
      do iw=1,nwave
        read(lu_rfi,'(a)') chars
        lw=lofw(chars)
        stations(iw) = chars(1:lw)
        write(fname(iw),'("rfi_files/OBS/",a,a3,a1)') chars(1:lw),
     1                                               '_ff',char(0)
        write(sname(iw),'("rfi_files/SYNT/",a,a3,a1)') chars(1:lw),
     1                                                '_ff',char(0)

        if(lroot)write(lu_out,*)
        if(lroot)write(lu_out,*) '* Read in data from ',
     1                            fname(iw)(1:lw+18)

        read(lu_rfi,*) weight(iw)
        if(lroot)write(lu_out,*) '* Weight =',weight(iw)
        if(lroot)write(lu_out,*)' '

        open(666, file=fname(iw), status="old")
        read(666,*) npoints
        read(666,*)
        do i=1,npoints
          read(666,*) freq_obs(iw,i),observed_data_r(i,iw),
     1                observed_data_i(i,iw)
          observed_data(i,iw) = complex(observed_data_r(i,iw),
     1                                  observed_data_i(i,iw))
        end do
        close(666)

        amp_obs_real(iw) = 0.0
        amp_obs_imag(iw) = 0.0
        do i=1, npoints
          amp_obs_real(iw) = amp_obs_real(iw) + real(observed_data(i,iw))**2
          amp_obs_imag(iw) = amp_obs_imag(iw) + imag(observed_data(i,iw))**2
        end do

      end do





      close(lu_rfi)


!       read kernels
      do i=1,nwave
        sta = stations(i)
        lw=lofw(sta)
!----------------------------------------------
!		load kernels first order
!----------------------------------------------
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

!----------------------------------------------
!		load kernels second order
!----------------------------------------------

        write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dS2dx2',char(0)
        open(11, file=filename, status="old")
        read(11, *) nlines(i)
        read(11,*)
        do j=1,nlines(i)
          read(11,*) dS2_freq(1,1,j), dS2_real(1,1,j), dS2_imag(1,1,j)
          dS2(i, 1,1,j)=complex(dS2_real(1,1,j),dS2_imag(1,1,j))

        end do
        close(11)

        write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdxdy',char(0)
        open(11, file=filename, status="old")
        read(11, *) nlines(i)
        read(11,*)
        do j=1,nlines(i)
          read(11,*) dS2_freq(1,2,j), dS2_real(1,2,j), dS2_imag(1,2,j)
          dS2(i, 1,2,j)=complex(dS2_real(1,2,j),dS2_imag(1,2,j))
        end do
        close(11)

        write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdxdz',char(0)
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
        open(11, file=filename, status="old")
        read(11, *) nlines(i)
        read(11,*)
        do j=1,nlines(i)
          read(11,*) dS2_freq(2,2,j), dS2_real(2,2,j), dS2_imag(2,2,j)
          dS2(i, 2,2,j)=complex(dS2_real(2,2,j),dS2_imag(2,2,j))
        end do
        close(11)

        write(filename, '("../src/rfi_subs/kernels/",a,a7,a1)') sta(1:lw),'_dSdydz',char(0)
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
        open(11, file=filename, status="old")
        read(11, *) nlines(i)
        read(11,*)
        do j=1,nlines(i)
          read(11,*) dS2_freq(3,3,j), dS2_real(3,3,j), dS2_imag(3,3,j)
          dS2(i, 3,3,j)=complex(dS2_real(3,3,j),dS2_imag(3,3,j))
        end do
        close(11)
!------------------------------------------------------------------------
!		read point source solution
!------------------------------------------------------------------------

        write(filename, '("../src/rfi_subs/point_source/",a,a3,a1)') sta(1:lw),'_ps',char(0)
        open(11, file=filename, status="old")
        read(11,*) nlines(i)
        read(11,*)
        do j=1, nlines(i)
          read(11,*) freq_ps(i, j), point_source_real(j), point_source_imag(j)
          point_source(i, j)=complex(point_source_real(j), point_source_imag(j))
        end do
        close(11)

      end do

      if(lroot)write(lu_out,*)

 100  format(1x,60("-")//)
 101  format(/1x,60("-")/)

      return
      end
!------------------------------------------------------------------------
      subroutine forward(nd,model,lppd)

      include 'rfi_subs/rfi_param.inc'

      real*4		lppd, misfitval

      real*4		model(nd)

      real*8		rmodel(max_nd)
      complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave)

      real*4		amp_obs_real(maxwave),
     &amp_obs_imag(maxwave)

      real*4		observed_data_r(maxdata,maxwave),
     &                  observed_data_i(maxdata,maxwave),
     &                  predicted_data_r(maxdata,maxwave),
     &                  predicted_data_i(maxdata,maxwave),
     &freq_obs(maxwave, maxdata),
     &weight(maxwave),
     &time_begin(maxwave),
     &time_end(maxwave),
     &Scalc_real(maxdata), Scalc_imag(maxdata)

      integer		ndata(maxwave), s
      integer		nwave, nd

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

!						Info and Logical unit common
!						blocks used by NA routines

      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &              lu_nad,verbose,debug,timing,summary


      common /rfi_com/observed_data, weight,
     &                  time_begin, time_end,  ndata,
     &                  nwave, lu_mod, fname, sname, stations,
     &dS1, dS2, point_source, freq_obs, nlines,
     &amp_obs_real, amp_obs_imag

!        Perform forward modelling on model `rmodel'.
      do j=1,nd
        rmodel(j) = dble(model(j))
      end do

      call forward_modelling(
     &    rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &dS1, dS2, point_source, freq_obs, nlines, .false., .false.)



!     Calculate misfit function

      call calcmisfit(
     &     predicted_data,
     &     observed_data, nlines,
     &     weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .false.)

      lppd =  misfitval * 1

      write(*,*) "Misfit: ", lppd
      write(123,*) "ooooooooooooooooooooo"
      write(123,*) "Tot Misfit: ", lppd
      write(123,*) "======================================"


      return
      end
!------------------------------------------------------------------------
      subroutine writemodels
     &             (nd, ntot, models, misfit, ns1, ns2, itmax,
     &              nh_max, nh, header)

      include 'rfi_subs/rfi_param.inc'

!						NA variables and arrays
      real*4		models(nd,*)
      real*4		misfit(ntot)
      real*4		mfitmin
      real*4		mfitminc
      real*4		mfitmean
      real*4		lppd, misfitval

      character*(*)   header

      logical		verbose,debug,timing,summary
      real*8		rmodel(max_nd)

      complex*8	observed_data(maxdata,maxwave),
     &                  predicted_data(maxdata,maxwave)

      real*4		amp_obs_real(maxwave),
     &amp_obs_imag(maxwave)

      real*4		weight(maxwave),
     &time_begin(maxwave),
     &time_end(maxwave),
     &freq_obs(maxwave, maxdata)

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

!						Info and Logical unit common
!						block used by NA routines

      common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,
     &              lu_nad,verbose,debug,timing,summary


      common /rfi_com/observed_data, weight,
     &                  time_begin, time_end,  ndata,
     &                  nwave, lu_mod, fname, sname, stations,
     &dS1, dS2, point_source, freq_obs, nlines,
     &amp_obs_real, amp_obs_imag


      nlayers = 6

!		write out models at each iteration

      mfitmin = misfit(1)
      ns = ns1
      np = 0
      mopt = 1

!		turn off writing to standard out by setting lu to zero

      lu_out2 = lu_out
      lu_out2 = 0
      write(lu_mod,*)ns1,' Number of samples in starting pool'
      write(lu_mod,*)ns2,' Number of new samples per iteration'
      write(lu_mod,*)itmax,' Number of iterations'
      write(*,*) "lu_mod: ", lu_mod

!		loop over iterations
      do it=1,itmax+1
        mfitminc = misfit(np+1)
        mfitmean = 0.0
!		find minimum and mean misfit
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

!			write from this iteration to output file.

        write(lu_mod,801) it-1, mfitmin, mfitmean, mfitminc
        do i=1,ns
          jj = np + i
          call display_model
     &         (lu_mod, i, models(1,jj), nd, misfit(jj) )
        end do
        np = np + ns
        ns = ns2

        call output_summary(
     &               lu_out2, lu_sum, it-1, models(1,mopt), nd,
     &               np, mfitmin, mfitmean, mfitminc, mopt)
      end do

!        Write out final model

      call display_final(
     &        lu_out, models(1,mopt), nd, mfitmin )

      call display_final(
     &        lu_sum, models(1,mopt), nd, mfitmin )

      call display_final(
     &        lu_mod, models(1,mopt), nd, mfitmin)

      do j=1,nd
        rmodel(j) = dble(models(j,mopt))
      end do

      write(*,*) "PEPPAAAAAAAAAAA"
      write(111,*) "************************  INITIAL MISFIT ************************ "
      call forward_modelling(
     &      rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &dS1, dS2, point_source, freq_obs, nlines, .false., .true.)


      call calcmisfit(
     &       predicted_data,
     &       observed_data, nlines,
     &       weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .true.)

      write(111,*) "***"
      write(111,*) "Total misfit: ", misfitval
      write(111,*) "***"

!	   repeat forward modelling for optimum model

      write(111,*) "************************ BEST MODEL ************************ "
      call forward_modelling(
     &    rmodel, nd, ndata, nwave, stations, predicted_data, observed_data,
     &dS1, dS2, point_source, freq_obs, nlines, .true., .false.)


      call calcmisfit(
     &     predicted_data,
     &     observed_data, nlines,
     &     weight, nwave, misfitval, fname, amp_obs_real, amp_obs_imag, .true.)

      write(111,*) "***"
      write(111,*) "total misfit: ", misfitval
      write(111,*) "***"
      close(111)

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

c	  write info into character string

      write(header(1:24),fmt='(4i6)')ns1,ns2,itmax,nlayers

      do i=1,nlayers
        k1 = 2*(i-1)*6 + 25
        k2 = k1 + 11
      end do


  801 format( 'iteration:',i5,',  misfit: min=',ES10.3,
     &        ', mean=',ES10.3,', minc=',ES10.3 )

      return
      end
