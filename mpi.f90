!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 module mpif
!*************************************************************************
   implicit none
   save
   
#ifdef _MPI
   include 'mpif.h'
   integer :: mpi_er, mpi_tg, mpi_rq(0:2*_Np), mpi_st(mpi_status_size)
#endif
   integer :: mpi_rnk, mpi_sze
   
 contains

!------------------------------------------------------------------------
!  initialise 
!------------------------------------------------------------------------
   subroutine mpi_precompute()
      mpi_rnk = 0
      mpi_sze = 1
#ifdef _MPI
      call mpi_init(mpi_er)
      call mpi_comm_rank(mpi_comm_world, mpi_rnk, mpi_er)
      call mpi_comm_size(mpi_comm_world, mpi_sze, mpi_er)
      if(mpi_sze /= _Np) stop 'mpi_precompute: incorrect num procs'
#endif
   end subroutine mpi_precompute

   
!*************************************************************************
 end module mpif
!*************************************************************************
