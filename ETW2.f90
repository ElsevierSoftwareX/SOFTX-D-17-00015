!*************************************************************************
! Energies as function of z
!  m=0  u_z     --> Emn
!  m=0  u_r,th  --> Esw
!  m/=0 u_z     --> Est
!  m/=0 u_r,th  --> Erl
!  omega_z      --> Ewz (enstrophy)
!*************************************************************************
#include "../parallel.h"
 PROGRAM MAIN
!*************************************************************************
   use io
   implicit none
   type (phys) :: p, p1,p2,p3
   double precision :: Emn(0:i_Z-1),Est(0:i_Z-1),Erl(0:i_Z-1),Esw(0:i_Z-1)
   double precision :: Ewz(0:i_Z-1),E__(0:i_Z-1)
   double precision :: zoffset, d1
   integer :: k, n
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
      print*, 'Enter state file:'
      read(*,9) io_statefile   
      print*, ' Enter shift data right zoffset:'
      read(*,*) zoffset
   else 
      open(99,status='old',file='ETW2.IN')
      read(99,9) io_statefile
      read(99,*) zoffset
      close(99)
   end if

   call io_load_state()
   call var_coll_shift(0d0,zoffset,vel_ur, vel_ur)
   call var_coll_shift(0d0,zoffset,vel_ut, vel_ut)
   call var_coll_shift(0d0,zoffset,vel_uz, vel_uz)
   call vel_transform()

!  m=0  u_z     --> Emn
!  m=0  u_r,th  --> Esw
!  m/=0 u_z     --> Est
!  m/=0 u_r,th  --> Erl
!  omega_z      --> Ewz (enstrophy)

   do n = 1, mes_D%pN
      do k = 0, i_Z-1
         p%Re(k,1,n) = sum(vel_r%Re(k,:,n))/dble(i_Th)
         p%Re(k,2,n) = sum(vel_t%Re(k,:,n))/dble(i_Th)
         p%Re(k,3,n) = sum(vel_z%Re(k,:,n))/dble(i_Th)
         p1%Re(k,:,n) = (vel_r%Re(k,:,n)-p%Re(k,1,n))**2
         p2%Re(k,:,n) = (vel_t%Re(k,:,n)-p%Re(k,2,n))**2
         p3%Re(k,:,n) = (vel_z%Re(k,:,n)-p%Re(k,3,n))**2

         p%Re(k,5,n) = p%Re(k,3,n)**2 * 2d0*d_PI
         p%Re(k,6,n) = (p%Re(k,1,n)**2+p%Re(k,2,n)**2) * 2d0*d_PI
         p%Re(k,7,n) = sum(p3%Re(k,:,n)) * 2d0*d_PI/dble(i_Th)
         p%Re(k,8,n) = sum(p1%Re(k,:,n)+p2%Re(k,:,n))*2d0*d_PI/dble(i_Th)

         p%Re(k,9,n) = sum(vel_curlz%Re(k,:,n)**2)*2d0*d_PI/dble(i_Th)
      end do   
   end do

   Emn = 0d0
   Esw = 0d0
   Est = 0d0
   Erl = 0d0
   Ewz = 0d0
   do n = 1, mes_D%pN
      do k = 0, i_Z-1
         Emn(k) = Emn(k) + 0.5d0*p%Re(k,5,n)*mes_D%intrdr(mes_D%pNi+n-1)
         Esw(k) = Esw(k) + 0.5d0*p%Re(k,6,n)*mes_D%intrdr(mes_D%pNi+n-1)
         Est(k) = Est(k) + 0.5d0*p%Re(k,7,n)*mes_D%intrdr(mes_D%pNi+n-1)
         Erl(k) = Erl(k) + 0.5d0*p%Re(k,8,n)*mes_D%intrdr(mes_D%pNi+n-1)
         Ewz(k) = Ewz(k) + 0.5d0*p%Re(k,9,n)*mes_D%intrdr(mes_D%pNi+n-1)
      end do
   end do
#ifdef _MPI      
   call mpi_allreduce(Emn, E__, i_Z, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   Emn = E__
   call mpi_allreduce(Esw, E__, i_Z, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   Esw = E__
   call mpi_allreduce(Est, E__, i_Z, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   Est = E__
   call mpi_allreduce(Erl, E__, i_Z, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   Erl = E__
   call mpi_allreduce(Ewz, E__, i_Z, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   Ewz = E__
#endif

   if(mpi_rnk==0) then
      print*, 'Writing Ez.dat'
      open(82, status='unknown', file='Ez.dat')   
      d1 = ((2d0*d_PI)/d_alpha) / dble(i_Z)
      write(82,*) '# Emn =',sum(Emn) * d1
      write(82,*) '# Esw =',sum(Esw) * d1
      write(82,*) '# Est =',sum(Est) * d1
      write(82,*) '# Erl =',sum(Erl) * d1
      write(82,*) '# Etot=',sum(Emn+Esw+Est+Erl) * d1
      do k = 0, i_Z-1
         d1 = ((2d0*d_PI)/d_alpha) * dble(k)/dble(i_Z)
         write(82,'(7e16.8)') d1, 0d0, Emn(k), Esw(k), Est(k), Erl(k), Ewz(k)
      end do
      close(82)
   end if
   
#ifdef _MPI
   call mpi_barrier(mpi_comm_world, mpi_er)
   call mpi_finalize(mpi_er)
#endif   

!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

