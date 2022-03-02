!------------------------------------------------------------------------
!
!    Subroutine read_vmodelNA - reads in parameters associated with
!                 velocity model and
!                 reads in all definitions to set up
!                 the parameter space
!
!    Note_1: a priori model co-variance (scale factors) for i-th
!        parameter given by scales(i+1). scales(1) determines
!        the type of relative axis scaling in parameter space
!        If scales(1) = 0 then; no scaling all other scale values ignored
!                     (This is only sensible when all parameters
!                  have the same physical units)
!        If scales(1) = -1 then; all parameter ranges are normalized
!                          and scales(j) j>1 are ignored
!                  (Effectively the parameter range becomes
!                   the a priori co-variance for each parameter)
!        If scales(1) = anything else; then scales(j+1) is taken as the
!                   a priori co-variance for parameter j.
!
!    Note_2: Calls no other routines
!
!       Note_3:
!               range(1,1:nlayer): lower bound of thickness
!               range(2,1:nlayer): upper bound of thickness
!               range(1,nlayer+1:2*nlayer): lower bound of velocity_1
!               range(2,nlayer+1:2*nlayer): upper bound of velocity_1
!               range(1,2*nlayer+1:3*nlayer): lower bound of velocity_2
!               range(2,2*nlayer+1:3*nlayer): upper bound of velocity_2
!                 velocity_1 : S-wave velocity at upper interface
!                 velocity_2 : S-wave velocity at lower interface
!               range(1,3*nlayer+1:4*nlayer): lower bound of Vp/Vs ratio
!               range(2,3*nlayer+1:4*nlayer): upper bound of Vp/Vs ratio
!
!-----------------------------------------------------------------------
!
      subroutine read_vmodelNA(
     1          lu_vel, range, scales, moddim)

         include 'homti_param.inc'

         real*4        range(2,*)
         real        scales(*)
         integer        moddim, nlayer

         read(lu_vel,*) nlayer,scales(1)
         write(*,*) "PAPA FRANCY", nlayer

         do i=1,nlayer
            i2=nlayer+i
            i3=2*nlayer+i
            read(lu_vel,*) range(1,i), range(2,i), scales(i+1),
     1                       range(1,i2), range(2,i2), scales(i2+1),
     1                       range(1,i3), range(2,i3), scales(i3+1)
         end do
         moddim=nlayer*3

         write(*,*) "SUOR MARIA CLARETTI", moddim

         return
      end

