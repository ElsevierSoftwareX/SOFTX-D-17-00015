#include "../parallel.h"
!*************************************************************************
 PROGRAM LOADSTATE
!*************************************************************************
   use io
   implicit none
   integer :: n
   _loop_km_vars
 9 format(A)

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

   print*, 'Enter input statefile:'
   read(*,9) io_statefile   
   
!  upload vel_ur, vel_ut, vel_uz, tim_t, tim_dt
!  interpolate onto new radial mesh if necessary
   call io_load_state()
   
!  get variables in physical space if required 
!    --> vel_r, vel_t, vel_z, vel_curlr, vel_curlt, vel_curlz
   call vel_transform()



   

!*************************************************************************
 END PROGRAM LOADSTATE
!*************************************************************************

