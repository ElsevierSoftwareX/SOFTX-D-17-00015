!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 module timestep
!*************************************************************************
   use variables
   implicit none
   save

   double precision :: tim_t
   double precision :: tim_dt
   double precision :: tim_dterr
   integer          :: tim_it
   integer          :: tim_step
   double precision :: tim_corr_dt
   double precision :: tim_cfl_dt
   integer          :: tim_cfl_dir
   logical          :: tim_new_dt
   integer          :: tim_pcn

 contains

!------------------------------------------------------------------------
!  set the initial time and timestep, 
!  possibly overwritten by loading a previously saved state
!------------------------------------------------------------------------
   subroutine tim_precompute()
      tim_t       = d_time
      tim_dt      = d_timestep
      tim_dterr   = 0d0
      tim_it      = 0
      tim_step    = 0
      tim_corr_dt = 0d0
      tim_cfl_dt  = 0d0
      tim_cfl_dir = 0
      tim_new_dt  = .false.
      tim_pcn     = -1
   end subroutine tim_precompute


!------------------------------------------------------------------------
!  Sets  A = c1 I + c2 Lap,   if PM=+-1 modified Lap for u+,u-
!------------------------------------------------------------------------
   subroutine tim_mesh_init(PM,c1,c2, A)
      integer,          intent(in)  :: PM
      double precision, intent(in)  :: c1,c2
      type (mesh),      intent(out) :: A(0:i_pH1)
      double precision :: d(i_N)
      _loop_km_vars
      
      _loop_km_begin
         d = -mes_D%r(:,-2)*i_Mp*m*i_Mp*m - d_alpha*k*d_alpha*k
         if(PM/=0) d = d - mes_D%r(:,-2) - 2d0*PM*i_Mp*m*mes_D%r(:,-2)
         A(nh)%M = c2 * mes_D%radLap%M
         A(nh)%M(i_KL+1,1:) = A(nh)%M(i_KL+1,1:) + c2*d + c1 
      _loop_km_end

   end subroutine tim_mesh_init


!------------------------------------------------------------------------
!  Sets  A = c1 I + c2 Lap(h)  then replaces A with its LU factorisation
!  if PM==0 usual laplacian, PM==-1/1 modified laplacian
!------------------------------------------------------------------------
   subroutine tim_lumesh_init(PM,BC,c1,c2, A)
      integer,          intent(in)  :: PM,BC
      double precision, intent(in)  :: c1,c2
      type (lumesh),    intent(out) :: A(0:i_pH1)
      double precision :: d(i_N)
      integer :: info, n,j, S
      _loop_km_vars
      
      _loop_km_begin
         d = -mes_D%r(:,-2)*i_Mp*m*i_Mp*m - d_alpha*k*d_alpha*k
         if(PM/=0) d = d - mes_D%r(:,-2) - 2d0*PM*i_Mp*m*mes_D%r(:,-2)
         A(nh)%M(i_KL+1:, :) = c2 * mes_D%radLap%M(:,1:)
         A(nh)%M(2*i_KL+1,:) = A(nh)%M(2*i_KL+1,:) + c2*d + c1 

         ! symmetry on axis            ! S==1 fn odd => f=0
         S = modulo(m*i_Mp+abs(PM),2)  ! S==0 fn even => dr(f)==0
         do j = 1, 1+i_KL
	    A(nh)%M(2*i_KL+1+1-j, j) = mes_D%dr0(j,1-S)
	 end do
					! boundary condition
         do j = i_N-i_KL, i_N
            A(nh)%M(2*i_KL+1+i_N-j,j) = mes_D%dr1(i_KL-i_N+j+1,BC) 
         end do
         if(BC==1 .and. k==0 .and. m==0) cycle

         call dgbtrf(i_N,i_N,i_KL,i_KL,A(nh)%M,3*i_KL+1,A(nh)%ipiv,info)
         if(info /= 0) stop 'tim_lumesh_init'
      _loop_km_end

   end subroutine tim_lumesh_init


!------------------------------------------------------------------------
!  solve system Ax=y for x, replaces b;  uses lapack routine dgbtrs()
!------------------------------------------------------------------------
   subroutine tim_lumesh_invert(BC,A, b)
      integer,       intent(in)    :: BC
      type (lumesh), intent(in)    :: A(0:i_pH1)
      type (coll),   intent(inout) :: b
      integer :: info
      _loop_km_vars

      _loop_km_begin
         if(BC==1 .and. k==0 .and. m==0) then
            b%Re(:,nh) = 0d0
            b%Im(:,nh) = 0d0         
            cycle
         end if 
         call dgbtrs('N', i_N, i_KL, i_KL, 1, A(nh)%M, 3*i_KL+1,  &
                     A(nh)%ipiv, b%Re(1,nh), i_N, info )
         if(info/=0) stop 'mes_lumesh_invert.1'
         call dgbtrs('N', i_N, i_KL, i_KL, 1, A(nh)%M, 3*i_KL+1,  &
                     A(nh)%ipiv, b%Im(1,nh), i_N, info )
         if(info/=0) stop 'mes_lumesh_invert.2'
      _loop_km_end
   
   end subroutine tim_lumesh_invert


!------------------------------------------------------------------------
!  multiply  d = A b + c  
!  S=0, b even for m even; S=1, b odd for m even
!------------------------------------------------------------------------
   subroutine tim_meshmult(S,A,b,c, d)
      integer,     intent(in)  :: S
      type (coll), intent(in)  :: b,c
      type (mesh), intent(in)  :: A(0:i_pH1)
      type (coll), intent(out) :: d
      double precision :: re(i_N), im(i_N)
      integer :: j,nl,nr,n
      _loop_km_vars

      _loop_km_begin
         re = c%Re(:,nh)
         im = c%Im(:,nh)
         do j = 1, i_N
            nl = max(1,j-i_KL)
            nr = min(j+i_KL,i_N)
            do n = nl, nr
               re(n) = re(n) + A(nh)%M(i_KL+1+n-j,j) * b%Re(j,nh)
               im(n) = im(n) + A(nh)%M(i_KL+1+n-j,j) * b%Im(j,nh)
            end do
         end do
         d%Re(:,nh) = re
         d%Im(:,nh) = im
      _loop_km_end
 
   end subroutine tim_meshmult


!-------------------------------------------------------------------------
!  set the RHS for the boundary condition = 0
!-------------------------------------------------------------------------
   subroutine tim_zerobc(a)
      type (coll), intent(out) :: a
      integer :: n
      do n = 0, var_H%pH1
         a%Re(  1, n ) = 0d0
         a%Re(i_N, n ) = 0d0
         a%Im(  1, n ) = 0d0
         a%Im(i_N, n ) = 0d0
      end do
   end subroutine tim_zerobc

!-------------------------------------------------------------------------
!  corrected "implicit" nonlinear terms
!-------------------------------------------------------------------------
   subroutine tim_nlincorr(N_, N)
      type (coll), intent(in)    :: N_
      type (coll), intent(inout) :: N
      double precision :: d1,d2
      integer :: nhm
      d1 = d_implicit
      d2 = 1d0 - d_implicit
      nhm = var_H%pH1
      N%Re(:,0:nhm) = d1 * N%Re(:,0:nhm) + d2 * N_%Re(:,0:nhm)
      N%Im(:,0:nhm) = d1 * N%Im(:,0:nhm) + d2 * N_%Im(:,0:nhm)
   end subroutine tim_nlincorr
      

!-------------------------------------------------------------------------
!  measure the magnutude of the correction: max diff
!-------------------------------------------------------------------------
   subroutine tim_measurecorr(c1,c2,c3)
      type (coll), intent(in)    :: c1,c2,c3
      double precision :: dterr, d
      integer :: nh
      dterr = 0d0
      do nh = 0, var_H%pH1
         dterr = max(dterr,  &
            maxval( c1%Re(:,nh)*c1%Re(:,nh) + c1%Im(:,nh)*c1%Im(:,nh)),  &
            maxval( c2%Re(:,nh)*c2%Re(:,nh) + c2%Im(:,nh)*c2%Im(:,nh)),  &
            maxval( c3%Re(:,nh)*c3%Re(:,nh) + c3%Im(:,nh)*c3%Im(:,nh)) )
      end do
#ifdef _MPI
      call mpi_allreduce(dterr, d, 1, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      dterr = d
#endif
      tim_dterr = max(tim_dterr,dsqrt(dterr))
   end subroutine tim_measurecorr


!-------------------------------------------------------------------------
!  check for convergence via the 2-norm of the correction
!  -- see tim_addcorr()
!-------------------------------------------------------------------------
   subroutine tim_check_cgce()
      double precision, save :: lasterr
      
      if(tim_it==1) then
         tim_corr_dt = tim_dt * dsqrt( d_dterr/tim_dterr )
         lasterr = 1d99
      end if         

      if(tim_dt<1d-9 .and. tim_step>30) then
         if(mpi_rnk==0) print*, 'tim_check_cgce: dt --> 0 !!!?'
         tim_step = i_maxtstep-1
         tim_it = 0
      else if(tim_it==tim_pcn) then
         tim_it = 0
         if(mpi_rnk==0 .and. modulo(tim_step,i_save_rate2)==0) then
            if(d_timestep> 0d0) print*,' step=',tim_step
            if(d_timestep<=0d0) print*,' step=',tim_step,' dt=',real(tim_dt)
         end if
      else if(tim_it>10) then
         if(mpi_rnk==0) print*, 'tim_check_cgce: too many its!!!'
         tim_step = i_maxtstep-1
         tim_it = 0
      else if(tim_dterr>lasterr) then
         if(mpi_rnk==0) print*, 'tim_check_cgce: increasing error!!!'
         if(tim_dterr>2d0*d_dterr) tim_step = i_maxtstep-1
         if(tim_dterr<2d0*d_dterr) tim_corr_dt = tim_dt/(1d1*d_courant)
         tim_it = 0
      else if(tim_dterr>d_dterr) then
         lasterr = tim_dterr
         tim_it = tim_it + 1
      else          
         if(mpi_rnk==0 .and. modulo(tim_step,i_save_rate2)==0) then
            if(d_timestep> 0d0) print*,' step=',tim_step,' its=',tim_it
            if(d_timestep<=0d0) print*,' step=',tim_step,' dt=',real(tim_dt)
         end if
         tim_it = 0
      end if
      tim_dterr = 0d0
      
   end subroutine tim_check_cgce


!-------------------------------------------------------------------------
!  set a new timestep size, if necessary.
!-------------------------------------------------------------------------
   subroutine tim_new_tstep()
      double precision :: dt
      integer, save :: i = 0

      dt = min(tim_dt*1.11d0, d_maxdt)
      if(tim_step==0d0 .and. tim_corr_dt==0d0)  &
                           dt = min(dt, tim_cfl_dt*0.1d0)
      if(tim_cfl_dt >0d0)  dt = min(dt, tim_cfl_dt*d_courant)
      if(tim_corr_dt>0d0)  dt = min(dt, tim_corr_dt*0.95d0)
      
      i = i - 1
      tim_new_dt = ( dt<tim_dt*0.95d0  &
              .or.  (dt>tim_dt*1.10d0 .and. i<0)  &
              .or.  (dt==d_maxdt .and. tim_dt<dt .and. i<0) )
      if(tim_new_dt) then
         if(dt>=tim_dt)  i = 99
         tim_dt = dt
      end if
      
   end subroutine tim_new_tstep


!*************************************************************************
 end module timestep
!*************************************************************************

