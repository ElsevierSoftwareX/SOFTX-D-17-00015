!*************************************************************************
#include "../parallel.h"
 PROGRAM ADDSTATES
!*************************************************************************
   use io
   implicit none
   double precision :: d1, E,Ek(0:i_K1),Em(0:i_M1)
   type (coll) :: r,t,z
   integer :: i, j
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
!   if(d_time>=0d0) stop 'set d_time<0 to upload t'
tim_t = 0d0
   
   print*, 'How many states?:  (Warning: writes state0000.cdf.dat)'
   read(*,*) i
   
   call var_coll_init(r)
   call var_coll_init(t)
   call var_coll_init(z)
   
   do j = 1, i
      print*, 'Enter state filename:'
      read(*,9) io_statefile
      call io_load_state()

      call var_coll_norm(vel_ur, d1,Ek,Em)
      E = d1
      call var_coll_norm(vel_ut, d1,Ek,Em)
      E = E + d1
      call var_coll_norm(vel_uz, d1,Ek,Em)
      E = E + d1
      print*, 'KE = ', E

      print*, 'Enter scale factor a :  (u_out += a * state)'
      read(*,*) d1
      r%Re = r%Re + d1*vel_ur%Re
      r%Im = r%Im + d1*vel_ur%Im
      t%Re = t%Re + d1*vel_ut%Re
      t%Im = t%Im + d1*vel_ut%Im
      z%Re = z%Re + d1*vel_uz%Re
      z%Im = z%Im + d1*vel_uz%Im
   end do

   call var_coll_copy(r, vel_ur)   
   call var_coll_copy(t, vel_ut)   
   call var_coll_copy(z, vel_uz)   
   
   call io_save_state()
   
!*************************************************************************
 END PROGRAM ADDSTATES
!*************************************************************************


