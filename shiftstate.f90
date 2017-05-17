!*************************************************************************
#include "../parallel.h"
 PROGRAM SHIFTSTATES
!*************************************************************************
   use io
   implicit none
   double precision :: th, z
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
   tim_t = 0d0
   
   print*, 'Enter state filename (Warning: writes state0000.cdf.dat):'
   read(*,9) io_statefile
   call io_load_state()

   print*, 'Enter shifts in theta and z :'
   read(*,*) th,z

   call var_coll_shift(th,z,vel_ur, vel_ur)
   call var_coll_shift(th,z,vel_ut, vel_ut)
   call var_coll_shift(th,z,vel_uz, vel_uz)
   
   call io_save_state()
   
!*************************************************************************
 END PROGRAM SHIFTSTATES
!*************************************************************************


