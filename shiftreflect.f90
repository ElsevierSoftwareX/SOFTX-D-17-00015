#include "goldensel.f90"
#include "../parallel.h"

 module rtz
   use variables
   save
   type (coll) :: r,t,z
   integer :: sym
 end module rtz


!*************************************************************************
 PROGRAM NORMDIFF
!*************************************************************************
   use rtz
   use io
   implicit none
   double precision :: a,u, f_,fmin,thmin,dth,th,fth
   external f_
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
   
   print*, 'Check which symmetry?'
   print*, '1, shift&reflect;  2, mirror'
   read(*,*) sym
   sym = 3 - sym     ! 1-->2, 2-->1; see vel_imposesym()

   print*, 'Enter state filename:'
   read(*,9) io_statefile
   call io_load_state()

   !  | u |
   call var_coll_copy(vel_ur, r)
   call var_coll_copy(vel_ut, t)
   call var_coll_copy(vel_uz, z)
   call norm(r,t,z, u)

   !  | u - u_s |
   fmin = 1d99
   th  = -d_PI*0.5d0
   dth = d_PI*0.01d0
   do while(th<d_PI*0.5d0)
      fth = f_(th)
      if(fth<fmin) thmin = th
      if(fth<fmin) fmin = fth
      th = th + dth   
   end do
   call goldensel(thmin-dth,thmin,thmin+dth,1d-15,f_, th,fth)
   
   print*, '       | u |         = ', real(u)
   print*, '   | u - u_s |      = ', real(fth)
   print*, '| u - u_s | / | u | = ', real(fth/u)
   print*, 'plane of symm at angle ', th
   print*
   print*, 'Write shifted state to state0000.cdf.dat?  1, true; 0 false'
   read(*,*) n
   if(n==1) then
      call var_coll_shift(th,0d0,vel_ur, vel_ur)
      call var_coll_shift(th,0d0,vel_ut, vel_ut)
      call var_coll_shift(th,0d0,vel_uz, vel_uz)
      call var_imposesym(sym, vel_ur,vel_ut,vel_uz)
      call io_save_state()
   end if
   

!*************************************************************************
 END PROGRAM NORMDIFF
!*************************************************************************

 subroutine norm(c1,c2,c3, E)
   use variables
   type (coll), intent(in) :: c1,c2,c3
   double precision, intent(out) :: E
   double precision :: E_, Ek(0:i_K1), Em(0:i_M1)

   call var_coll_norm(c1, E ,Ek,Em)
   call var_coll_norm(c2, E_,Ek,Em)
   E   = E + E_
   call var_coll_norm(c3, E_,Ek,Em)
   E   = E + E_
   E = dsqrt(2d0*E)
      
 end subroutine norm


 double precision function f_(th)   
   use rtz
   implicit none
   double precision, intent(in)  :: th
   type (coll) :: c1,c2,c3, c4,c5,c6
   double precision :: a
   call var_coll_shift(th,0d0,r, c1)
   call var_coll_shift(th,0d0,t, c2)
   call var_coll_shift(th,0d0,z, c3)
   call var_coll_copy(c1, c4)
   call var_coll_copy(c2, c5)
   call var_coll_copy(c3, c6)
   call var_imposesym(sym, c4,c5,c6)
   call var_coll_sub(c4, c1)
   call var_coll_sub(c5, c2)
   call var_coll_sub(c6, c3)
   call norm(c1,c2,c3, a)
   f_ = a
 end function f_
 

