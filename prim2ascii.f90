#include "../parallel.h"
!*************************************************************************
 PROGRAM MAIN
!*************************************************************************
   use io
   implicit none
   double precision :: dRe, alpha, c
   integer :: e,f,i, n
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
   
   if(mpi_sze/=1) stop 'set _Np=1'
   
   print*, 'Enter input statefile:'
   read(*,9) io_statefile   

   print*, 'Include phase speed in output?:  0 false, 1 true'
   read(*,*) i
   if(i==1) then
      print*, ' enter c:'
      read(*,*) c
   end if

   print*, 'loading state...'
   call io_load_state()
   vel_uz%Re(:,0) = vel_uz%Re(:,0) + vel_U

   e=nf90_open(io_statefile,nf90_nowrite, f)
   e=nf90_get_att(f,nf90_global,'Re',    dRe)
   e=nf90_get_att(f,nf90_global,'alpha', alpha)
   e=nf90_close(f)
   
   print*, 'writing state.dat...'
   open(10,status='unknown',file='state.dat')
   if(i==1) write(10,'(2e20.12,I3,1e20.12)') dRe, alpha, i_Mp, c
   if(i/=1) write(10,'(2e20.12,I3)') dRe, alpha, i_Mp
   write(10,*) i_N, i_K, i_M
   do n = 1, i_N
      write(10,'(1e16.8)') mes_D%r(n,1)
   end do
   _loop_km_begin
     do n = 1, i_N
        write(10,'(2I4,6e16.8)')  &
           k, m, &
           vel_ur%Re(n,nh), vel_ur%Im(n,nh),  &
           vel_ut%Re(n,nh), vel_ut%Im(n,nh),  &
           vel_uz%Re(n,nh), vel_uz%Im(n,nh)
     end do
   _loop_km_end
   close(10)

        
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

