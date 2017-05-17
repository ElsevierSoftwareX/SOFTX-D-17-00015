#include "../parallel.h"
!*************************************************************************
 PROGRAM LOADSTATE
!*************************************************************************
   use io
   implicit none
   double precision :: mean_u(i_N), stdv_u(i_N)
   double precision :: mean_v(i_N), stdv_v(i_N)
   double precision :: mean_w(i_N), stdv_w(i_N)
   double precision :: d1, d(i_N)
   character(4) :: cnum
   integer :: n,n_, ff,fl,f
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
   
   if(mpi_sze==1) then
      print*, 'Enter first and last file numbers'
      print*, ' (current dir; approx equally spaced in time):'
      read(*,*)  ff, fl
   else
      open(99,status='old',file='MEANFLUCT.IN')
      read(99,*) ff, fl
      close(99)
   end if

   mean_u = 0d0
   stdv_u = 0d0
   mean_v = 0d0
   stdv_v = 0d0
   mean_w = 0d0
   stdv_w = 0d0

   do f = ff, fl
      write(cnum,'(I4.4)') f
      io_statefile = 'state'//cnum//'.cdf.dat'
      call io_load_state()
      call vel_transform()
      do n = 1, mes_D%pN
         n_ = mes_D%pNi + n - 1
         mean_u(n_) = mean_u(n_) + sum(vel_r%Re(:,:,n))
         stdv_u(n_) = stdv_u(n_) + sum(vel_r%Re(:,:,n)**2)
         mean_v(n_) = mean_v(n_) + sum(vel_t%Re(:,:,n))
         stdv_v(n_) = stdv_v(n_) + sum(vel_t%Re(:,:,n)**2)
         mean_w(n_) = mean_w(n_) + sum(vel_z%Re(:,:,n))
         stdv_w(n_) = stdv_w(n_) + sum(vel_z%Re(:,:,n)**2)
      end do
   end do

#ifdef _MPI
   call mpi_allreduce(mean_u, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_u = d
   call mpi_allreduce(stdv_u, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_u = d
   call mpi_allreduce(mean_v, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_v = d
   call mpi_allreduce(stdv_v, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_v = d
   call mpi_allreduce(mean_w, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   mean_w = d
   call mpi_allreduce(stdv_w, d, i_N, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   stdv_w = d
#endif

   d1 = dble((fl-ff+1) * i_Th*i_Z)
   mean_u = mean_u / d1
   stdv_u = dsqrt( stdv_u/d1 - mean_u**2 )
   mean_v = mean_v / d1
   stdv_v = dsqrt( stdv_v/d1 - mean_v**2 )
   mean_w = mean_w / d1
   stdv_w = dsqrt( stdv_w/d1 - mean_w**2 )
   mean_w = mean_w + 1d0-mes_D%r(:,2)
   
   if(mpi_rnk==0) then
      print*, 'Writing vel_meanstdv.dat ...'
      open(99, status='unknown', file='vel_meanstdv.dat')
      do n = 1, i_N
         write(99,'(7e16.8)') mes_D%r(n,1),  &
            mean_u(n), stdv_u(n), mean_v(n), stdv_v(n), mean_w(n), stdv_w(n)
      end do
      close(99)
   endif

#ifdef _MPI
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)
#endif

!*************************************************************************
 END PROGRAM LOADSTATE
!*************************************************************************

