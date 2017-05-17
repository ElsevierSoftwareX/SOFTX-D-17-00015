#include "../parallel.h"
!*************************************************************************
 PROGRAM IC
!*************************************************************************
   use io
   implicit none
   type (spec) :: s
   double precision :: r,th,z
   integer :: n
   _loop_km_vars

   print*, 'initialising...'
   call mpi_precompute()
   call par_precompute()
   call mes_precompute()
   call var_precompute()
   call tra_precompute()
   call tim_precompute()
   call vel_precompute()
   call  io_precompute()
   
   if(mpi_sze/=1) stop 'set _Np 1'

   tim_t  = 0d0
   tim_dt = 0.001d0
   
! **** SET IC IN REAL SPACE ****
   do n = 1, i_N
      do m = 0, i_Th-1
         do k = 0, i_Z-1
            r  = mes_D%r(n,1)
            th = (2d0*d_PI/dble(i_Th)) * m
            z  = (2d0*d_PI/(d_alpha*i_Z)) * k
            vel_r%Re(k,m,n) =    
            vel_t%Re(k,m,n) =   
            vel_z%Re(k,m,n) =   
         end do   
      end do
   end do
   call tra_phys2spec(vel_r, s)
   call var_spec2coll(s, vel_ur)
   call tra_phys2spec(vel_t, s)
   call var_spec2coll(s, vel_ut)
   call tra_phys2spec(vel_z, s)
   call var_spec2coll(s, vel_uz)

! **** SET IC IN COLLOCATED SPACE ****
! k and m are set by the _loop_km macro,
! see parallel.h

   _loop_km_begin
      do n = 1, i_N
         r = mes_D%r(n,1)
         vel_ur%Re(n,nh) =   
         vel_ur%Im(n,nh) =   
         vel_ut%Re(n,nh) =   
         vel_ut%Im(n,nh) =   
         vel_uz%Re(n,nh) =   
         vel_uz%Im(n,nh) =   
      end do
   _loop_km_end

! **** SAVE TO state0000.cdf.dat ****
   call io_save_state()


!*************************************************************************
 END PROGRAM IC
!*************************************************************************

