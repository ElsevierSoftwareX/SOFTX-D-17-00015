!*************************************************************************
#ifndef main_f90
#define main_f90
#define var_null util
#include "../program/main.f90"
#endif
!*************************************************************************
 module nonnewtonian
   use variables
   implicit none
   save
   
   double precision :: rho, radius		! (kg/m^3), (m)

   integer          :: model
						! Carreau-Yasuda model
   double precision :: lam_CY, n_CY, a_CY	! lam (s)
   double precision :: mu_0, mu_inf		! (Pa s = kg/(m.s))
   
   double precision :: fc1			! coeff of applied force
   
   type (coll) :: ur, ut, uz, nuc
   type (phys) :: S(3,3), sr, nu
   double precision :: nu_wall, nu_min, nu_max
   double precision :: lam_			! non-dim lam
   
 end module nonnewtonian   

!-------------------------------------------------------------------------
! called from MAIN
!-------------------------------------------------------------------------
 subroutine util(FLAG)
   use io
   implicit none
   integer, intent(in) :: FLAG
   logical, save :: set = .false.
   
   if(.not.set) then
      set = .true.
      call init_params_model()
   end if
				! called just before writing data in main   
   if(FLAG==1) then
      if(modulo(tim_step,i_save_rate1)==0) then
         call save_profiles()
      end if
      if(modulo(tim_step,i_save_rate2)==0) then
         call write_Re_nu()
      end if
      				! called after eval of nonlin terms
   else if(FLAG==2) then
      call get_nu()
      call update_nonlin_terms()
      
   end if
      
 end subroutine util

!-------------------------------------------------------------------------
! load params
!-------------------------------------------------------------------------
 subroutine init_params_model()
   use nonnewtonian
   use timestep
   implicit none

   tim_pcn = 1
   tim_it = 0

   open(99,status='old',file='params_model.in')
   read(99,*)  rho, radius
   read(99,*)  model
   if(model==1) then
      if(mpi_rnk==0) print*, 'Loading Carreau-Yasuda parameters...'
      read(99,*) mu_0, mu_inf
      read(99,*) lam_CY, n_CY, a_CY
      lam_ = d_Re * lam_CY * mu_0 / (rho * radius**2)
   else
      stop 'load_params_model: unknown model' 
   end if
   read(99,*)  fc1
   close(99)

   call get_nu()

 end subroutine init_params_model

   
!-------------------------------------------------------------------------
! calc strain rate, nu, nu_max/min/wall
!-------------------------------------------------------------------------
 subroutine get_nu()
   use nonnewtonian
   use velocity
   implicit none
   double precision :: d1,d2
   type (coll) :: c1,c2,c3
   type (spec) :: s1,s2,s3
   integer :: m_,n, i, j
   _loop_km_vars

				   ! vec{u} = (ur,ut,uz)
   call var_coll_copy(vel_ur, ur)
   call var_coll_copy(vel_ut, ut)
   call var_coll_copy(vel_uz, uz)
   if(mpi_rnk==0) uz%Re(:,0) = uz%Re(:,0) + vel_U

				   ! get velocity gradient tensor A
   call var_coll_meshmult(1,mes_D%dr(1),ur, c1)
   call var_coll_meshmult(1,mes_D%dr(1),ut, c2)
   call var_coll_meshmult(0,mes_D%dr(1),uz, c3)
   call var_coll2spec(c1,s1, c2,s2, c3,s3)
   call tra_spec2phys(s1, S(1,1))
   call tra_spec2phys(s2, S(1,2))
   call tra_spec2phys(s3, S(1,3))
   _loop_km_begin
      m_ = m*i_Mp
      c1%Re(:,nh) = mes_D%r(:,-1)*(-ur%Im(:,nh)*m_-ut%Re(:,nh))
      c1%Im(:,nh) = mes_D%r(:,-1)*( ur%Re(:,nh)*m_-ut%Im(:,nh))
      c2%Re(:,nh) = mes_D%r(:,-1)*(-ut%Im(:,nh)*m_+ur%Re(:,nh))
      c2%Im(:,nh) = mes_D%r(:,-1)*( ut%Re(:,nh)*m_+ur%Im(:,nh))
      c3%Re(:,nh) = mes_D%r(:,-1)*(-uz%Im(:,nh)*m_)
      c3%Im(:,nh) = mes_D%r(:,-1)*( uz%Re(:,nh)*m_)
   _loop_km_end
   call var_coll2spec(c1,s1, c2,s2, c3,s3)
   call tra_spec2phys(s1, S(2,1))
   call tra_spec2phys(s2, S(2,2))
   call tra_spec2phys(s3, S(2,3))
   _loop_km_begin
      d1 = d_alpha*k
      c1%Re(:,nh) = -ur%Im(:,nh)*d1
      c1%Im(:,nh) =  ur%Re(:,nh)*d1
      c2%Re(:,nh) = -ut%Im(:,nh)*d1
      c2%Im(:,nh) =  ut%Re(:,nh)*d1
      c3%Re(:,nh) = -uz%Im(:,nh)*d1
      c3%Im(:,nh) =  uz%Re(:,nh)*d1
   _loop_km_end
   call var_coll2spec(c1,s1, c2,s2, c3,s3)
   call tra_spec2phys(s1, S(3,1))
   call tra_spec2phys(s2, S(3,2))
   call tra_spec2phys(s3, S(3,3))
				   ! get S = (1/2)(A + A^T)
   S(1,2)%Re = 0.5d0 *(S(1,2)%Re+S(2,1)%Re)
   S(1,3)%Re = 0.5d0 *(S(1,3)%Re+S(3,1)%Re)
   S(2,3)%Re = 0.5d0 *(S(2,3)%Re+S(3,2)%Re)
   S(2,1)%Re = S(1,2)%Re
   S(3,1)%Re = S(1,3)%Re
   S(3,2)%Re = S(2,3)%Re
				   ! strain rate = (2 S_ij S_ij)^(1/2)
   sr%Re =    S(1,1)%Re**2 + S(2,2)%Re**2 + S(3,3)%Re**2  &
      + 2d0*( S(1,2)%Re**2 + S(1,3)%Re**2 + S(2,3)%Re**2 )
   sr%Re = dsqrt( 2d0 * sr%Re )
				   ! get nu   
   if(model==1) then  ! Carreau-Yasuda
      d1 = mu_inf/mu_0
      d2 = - n_CY/a_CY
      nu%Re = d1 + (1d0-d1)*(1d0+(lam_*sr%Re)**a_CY)**d2
   else 
      stop 'get_nu: model not implemented'
   end if
   call tra_phys2spec(nu, s1)
   call var_spec2coll(s1, nuc)
   				   ! nu_wall,min,max
   nu_wall = nuc%Re(i_N,0)
   nu_min = minval(nu%Re)
   nu_max = maxval(nu%Re)
#ifdef _MPI
   call mpi_bcast(nu_wall, 1, mpi_double_precision,  &
      0, mpi_comm_world, mpi_er)
   call mpi_allreduce(nu_min, d1, 1, mpi_double_precision,  &
      mpi_min, mpi_comm_world, mpi_er)
   nu_min = d1
   call mpi_allreduce(nu_max, d1, 1, mpi_double_precision,  &
      mpi_max, mpi_comm_world, mpi_er)
   nu_max = d1
#endif

   if(tim_step/=0) return
   vel_nu = 1d0

/*
   ! Reference viscosity vel_nu should be approx equal to largest
   ! observed viscosity.  If it deviates to far, update vel_nu
   ! and recalc vel_matrices().
   ! For the shear-thinning model, we use nu_max in the Reynolds number
   ! so that vel_nu=1 and does not change.
   if(tim_it/=0) return
   n_ = mes_D%pN
   nu_max = maxval(nu%Re(:,:,1:n_))
#ifdef _MPI
   call mpi_allreduce(nu_max, d1, 1, mpi_double_precision,  &
      mpi_max, mpi_comm_world, mpi_er)
   nu_max = d1
#endif
   if(0.95d0*d1<vel_nu .and. vel_nu<1.05d0*d1)  return
   vel_nu = d1
   call vel_matrices()
*/
 end subroutine get_nu

!-------------------------------------------------------------------------
! add in terms from variable viscosity
!-------------------------------------------------------------------------
 subroutine update_nonlin_terms()
   use nonnewtonian
   use velocity
   implicit none
   type (coll) :: c1,c2,c3
   type (spec) :: s1,s2,s3
   type (phys) :: q(3), g(3), nu_
   logical, save :: set = .false.
   integer :: n, i, j
   _loop_km_vars
						   ! grad(nu')
   call var_coll_grad(nuc, c1,c2,c3)
   call var_coll2spec(c1,s1, c2,s2, c3,s3)
   call tra_spec2phys(s1, g(1))
   call tra_spec2phys(s2, g(2))
   call tra_spec2phys(s3, g(3))
						   ! curl(curl(u))
   call var_coll_curl(ur,ut,uz, c1,c2,c3)
   call var_coll_curl(c1,c2,c3, c1,c2,c3)
   call var_coll2spec(c1,s1, c2,s2, c3,s3)
   call tra_spec2phys(s1, q(1))
   call tra_spec2phys(s2, q(2))
   call tra_spec2phys(s3, q(3))

   ! N' = (1/Re)[(nu-nu_r)Lap(u) + grad(nu).2S + (nu_r-1)(-4\hat{z})]
   nu_%Re = nu%Re - vel_nu
   do i = 1, 3
      q(i)%Re = nu_%Re*(-q(i)%Re)
      do j = 1, 3
         q(i)%Re = q(i)%Re  +  g(j)%Re * S(j,i)%Re * 2d0
      end do
      if(i==3) q(i)%Re = q(i)%Re + (vel_nu-1d0)*(-4d0) 
      q(i)%Re = q(i)%Re / d_Re
   end do      
   call tra_phys2spec(q(1), s1)
   call tra_phys2spec(q(2), s2)
   call tra_phys2spec(q(3), s3)
   call var_spec2coll(s1,c1, s2,c2, s3,c3)
						   ! N += N'         
   call var_coll_add(c1, vel_Nr)
   call var_coll_add(c2, vel_Nt)
   call var_coll_add(c3, vel_Nz)
   
 end subroutine update_nonlin_terms


!-------------------------------------------------------------------------
! write to files for non-newtonian data
!-------------------------------------------------------------------------
 subroutine save_profiles()
   use nonnewtonian
   use io
   implicit none
   double precision :: r(i_N), vz(-i_N:i_N), vz_(-i_N:i_N)
   character(4) :: cnum
   integer :: i,n, m1,m2,m3, pN,N1,N2

   r = mes_D%r(:,1)
   write(cnum,'(I4.4)') io_save1

   if(mpi_rnk==0) then
      open(11, status='unknown', file='nu_prof'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# r  nu(r)'
      write(11,'(2e20.12)') (r(n), nuc%Re(n,0), n=1,i_N)
      close(11)
   end if

   if(_Ns/=1) stop 'write2files: err1'
   m1 = i_Th/4
   m2 = (2*i_Th)/4
   m3 = (3*i_Th)/4
   pN = mes_D%pN
   N1 = mes_D%pNi
   N2 = mes_D%pNi + pN - 1
   do i = 1, 2
      vz = -1d99
      if(i==1) vz(-N1:-N2:-1) = vel_z%Re(0,m2,1:pN) + vel_U(N1:N2)
      if(i==1) vz( N1: N2)    = vel_z%Re(0, 0,1:pN) + vel_U(N1:N2)
      if(i==2) vz(-N1:-N2:-1) = vel_z%Re(0,m3,1:pN) + vel_U(N1:N2)
      if(i==2) vz( N1: N2)    = vel_z%Re(0,m1,1:pN) + vel_U(N1:N2)
#ifdef _MPI
      call mpi_allreduce(vz, vz_, 2*i_N+1, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      vz = vz_
#endif
      if(mpi_rnk==0) then
        if(i==1) open(11,status='unknown',file='vel_prof'//cnum//'_0.dat')
        if(i==2) open(11,status='unknown',file='vel_prof'//cnum//'_1.dat')
        write(11,*) '# t = ', tim_t
        write(11,*) '# r  uz(r)'
        write(11,'(2e20.12)') (-r(-n), vz(n), n=-i_N,-1)
        write(11,'(2e20.12)') ( r( n), vz(n), n=1,i_N)
        close(11)
      end if
   end do

 end subroutine save_profiles


!-------------------------------------------------------------------------
! write to files for non-newtonian data
!-------------------------------------------------------------------------
 subroutine write_Re_nu()
   use nonnewtonian
   use io
   implicit none
   logical, save :: set=.false.

   if(mpi_rnk/=0) return

   if(.not.set) then
      set = .true.
      open(54,status='unknown',file='nu_minmax.dat')
   end if
   
   write(54,'(5e16.8)') tim_t, vel_nu, nu_wall, nu_min, nu_max
   
 end subroutine write_Re_nu
