!-------------------------------------------------------------------------
!
!    Subroutine summary - writes out a summary of the current
!                 iteration of the Neighbourhood algorithm
!                 together with the model `rmodel'.
!
!      Modified from original NA package
!--------------------------------------------------------------------------
!
      subroutine output_summary(
     1lu_out, lu_sum, it, rmodel, moddim, ntot,
     1mfitmin, mfitmean, mfitminc, mopt)

         include 'homti_param.inc'


         real*4          rmodel(*)

         real*4          mfitmin,
     1                   mfitmean,
     1                   mfitminc
         logical         lw

!                       Write out headers for files
!                       with summary of optimization
!                       performance.
         lw = .true.
         if(lu_out.eq.0)lw = .false.
         write(lu_sum,801)
         if(lw)write(lu_out,801)

         write(lu_sum,811) it,ntot,mfitmin,mfitmean,
     1                       mfitminc,mopt
         if(lw)write(lu_out,811) it,ntot,mfitmin,mfitmean,
     1                       mfitminc,mopt

         nlayer=moddim/4
         do j=1,nlayer
            j2=nlayer+j
            j3=2*nlayer+j
            j4=3*nlayer+j
            write(lu_sum,812) j, rmodel(j), rmodel(j2),
     1                          rmodel(j3), rmodel(j4)
            if(lw)write(lu_out,812) j, rmodel(j), rmodel(j2),
     1                          rmodel(j3), rmodel(j4)
         end do

         write(lu_sum,*)
         if(lw)write(lu_out,*)

  801    format( 3x,'It',1x,'Nsampled',3x,'Mfitmin',
     1             2x,'Mfitmean',1x,'Mfitmeanc',1x,
     1                3x,'Mopt')
  811    format( i5,i9,3f10.5,1x,i7)
  812    format( 5x,i3,4f10.3,1x,i7)

         return
      end
