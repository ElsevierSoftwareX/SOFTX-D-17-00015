!*************************************************************************
! Typical params_model.in:
!
! Smagorinsky:
!     1
!   0.025  1  25.      
!
! Vreman:
!     2
!   0.03   2  0.
!
!*************************************************************************
#ifndef main_f90
#define main_f90
#define var_null util
#include "../program/main.f90"
#endif
#include "bisect.f90"
!*************************************************************************
 module les
   use variables
   implicit none
   save
   
   integer          :: model, filter
   double precision :: C_s
   double precision :: A_plus
   
   double precision :: Re_tau, u_tau, nu_max
   double precision :: Rdth, dz, Delta_r(i_N), y_plus(i_N), vanDr(i_N)
   type (coll) :: ur, ut, uz, nuc
   type (phys) :: A(3,3), S(3,3), sr, nu
   
 end module les

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
      call bisect(FLAG)
      if(modulo(tim_step,i_save_rate1)==0)  &
         call save_nuprof()
      if(modulo(tim_step,i_save_rate2)==0)  &
         call write_Retau()
      				! called after eval of nonlin terms
   else if(FLAG==2) then
      call get_nu()
      if(tim_it==0) call update_nu_r()
      call update_nonlin_terms()
      
   end if
      
 end subroutine util

!-------------------------------------------------------------------------
! load params
!-------------------------------------------------------------------------
 subroutine init_params_model()
   use les
   use timestep
   implicit none
   integer :: n
 
   tim_pcn = 1
   tim_it  = 0

   open(99,status='old',file='params_model.in')
   read(99,*)  model
   read(99,*)  C_s, filter, A_plus
   close(99)

   Rdth = (2d0*d_PI/dble(i_Mp*i_Th)) * (3d0/2d0)  
   dz   = (2d0*d_PI/(d_alpha*i_Z))   * (3d0/2d0)
   Delta_r(1) = (mes_D%r(2,1)-0d0)/2d0
   do n = 2, i_N-1
      Delta_r(n) = (mes_D%r(n+1,1)-mes_D%r(n-1,1))/2d0
   end do
   Delta_r(i_N) = 1d0-mes_D%r(i_N-1,1)
   if(mpi_rnk==0) then
      print*, ' LES: max(dr) = ', maxval(Delta_r)
      print*, ' LES:   R.dth = ', Rdth
      print*, ' LES:      dz = ', dz
   end if      

 end subroutine init_params_model

   
!-------------------------------------------------------------------------
! calc strain rate, nu, nu_max/min/wall
!-------------------------------------------------------------------------
 subroutine get_nu()
   use les
   use velocity
   implicit none
   double precision :: d1,d2
   type (coll) :: c1,c2,c3
   type (spec) :: s1,s2,s3
   integer :: m_,n,n_, i,j
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
   call tra_spec2phys(s1, A(1,1))
   call tra_spec2phys(s2, A(1,2))
   call tra_spec2phys(s3, A(1,3))
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
   call tra_spec2phys(s1, A(2,1))
   call tra_spec2phys(s2, A(2,2))
   call tra_spec2phys(s3, A(2,3))
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
   call tra_spec2phys(s1, A(3,1))
   call tra_spec2phys(s2, A(3,2))
   call tra_spec2phys(s3, A(3,3))
				   ! get S = (1/2)(A + A^T)
   S(1,1)%Re = A(1,1)%Re
   S(2,2)%Re = A(2,2)%Re
   S(3,3)%Re = A(3,3)%Re
   S(1,2)%Re = 0.5d0 *(A(1,2)%Re+A(2,1)%Re)
   S(1,3)%Re = 0.5d0 *(A(1,3)%Re+A(3,1)%Re)
   S(2,3)%Re = 0.5d0 *(A(2,3)%Re+A(3,2)%Re)
   S(2,1)%Re = S(1,2)%Re
   S(3,1)%Re = S(1,3)%Re
   S(3,2)%Re = S(2,3)%Re
				   ! strain rate = (2 S_ij S_ij)^(1/2)
   sr%Re =    S(1,1)%Re**2 + S(2,2)%Re**2 + S(3,3)%Re**2  &
      + 2d0*( S(1,2)%Re**2 + S(1,3)%Re**2 + S(2,3)%Re**2 )
   sr%Re = dsqrt( 2d0 * sr%Re )
				   ! u_tau, Re_tau, y^+
   d1    = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
   u_tau = dsqrt(dabs((d1-2d0)/d_Re))
#ifdef _MPI
   call mpi_bcast(u_tau, 1, mpi_double_precision,  &
      0, mpi_comm_world, mpi_er)
#endif
   Re_tau = d_Re * u_tau
   y_plus = Re_tau*(1d0-mes_D%r(:,1))
				   ! eddy visc model
   if(model==1) call Smagorinsky()
   if(model==2) call Vreman()
   call tra_phys2spec(nu, s1)
   call var_spec2coll(s1, nuc)
				   ! nu_max
   n_ = mes_D%pN
   nu_max = maxval(nu%Re(:,:,1:n_))
#ifdef _MPI
   call mpi_allreduce(nu_max, d1, 1, mpi_double_precision,  &
      mpi_max, mpi_comm_world, mpi_er)
   nu_max = d1
#endif

 end subroutine get_nu


!-------------------------------------------------------------------------
! eddy viscosity model
!-------------------------------------------------------------------------
 subroutine Smagorinsky()
   use les
   use velocity
   implicit none
   double precision :: d1
   integer :: n, n_
   double precision, save :: Delta(i_N)
   logical, save :: set = .false.

   if(.not.set) then   
      set = .true.
      if(filter==1) then
         Delta = (Delta_r * Rdth * dz)**(1d0/3d0)
      else if(filter==2) then
         Delta = (Delta_r * mes_D%r(:,1)*Rdth * dz)**(1d0/3d0)
      else if(filter==3) then
         Delta = sqrt((Delta_r**2 + (mes_D%r(:,1)*Rdth)**2 + dz**2)/3d0)
      end if
   end if
   				   ! vanDriest fn
   vanDr  = 1d0 - dexp(-(y_plus/A_plus)**2)
				   ! nu_t
   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      d1 = vanDr(n_) * (C_s*Delta(n_))**2
      nu%Re(:,:,n) = d1 * sr%Re(:,:,n)
   end do
   
 end subroutine Smagorinsky


!-------------------------------------------------------------------------
! eddy viscosity model
!-------------------------------------------------------------------------
 subroutine Vreman()
   use les
   use velocity
   implicit none
   type (phys) :: B(3,3), P, AA
   double precision :: d1
   integer :: i,j,m, n,n_

   AA%Re = 0d0
   do i = 1, 3
      do j = 1, 3
         AA%Re = AA%Re + A(i,j)%Re**2
         B(i,j)%Re = 0d0
         do m = 1, 3
            do n = 1, mes_D%pN
               n_ = mes_D%pNi + n - 1
               if(m==1) d1 = Delta_r(n_)
               if(m==2 .and. filter==1) d1 = Rdth
               if(m==2 .and. filter==2) d1 = mes_D%r(n_,1)*Rdth
               if(m==3) d1 = dz
               B(i,j)%Re(:,:,n) = B(i,j)%Re(:,:,n)  &
                  + (d1**2) * A(m,i)%Re(:,:,n) * A(m,j)%Re(:,:,n)
            end do         
         end do
      end do
   end do
   
   P%Re = B(1,1)%Re*B(2,2)%Re - B(1,2)%Re**2  &
        + B(1,1)%Re*B(3,3)%Re - B(1,3)%Re**2  &
        + B(2,2)%Re*B(3,3)%Re - B(2,3)%Re**2
   P%Re  = merge(0d0,P%Re,AA%Re<1d-12)
   P%Re  = max(P%Re,0d0)
   AA%Re = max(AA%Re,1d-12)
   
   nu%Re = C_s * dsqrt(P%Re/AA%Re)

   if(A_plus<=0d0) return
   vanDr  = 1d0 - dexp(-(y_plus/A_plus)**2)
   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      d1 = 0.03d0/C_s
      d1 = d1 + (1d0-d1)* vanDr(n_)
      nu%Re(:,:,n) = d1 * nu%Re(:,:,n)
   end do
 
 end subroutine Vreman


!-------------------------------------------------------------------------
! check if the reference viscosity needs changing
!-------------------------------------------------------------------------
 subroutine update_nu_r()
   use les
   use velocity
   double precision :: d1
   d1 = 1d0 + d_Re*nu_max
   d1 = 1.10d0*d1
   if(vel_nu<0.93d0*d1 .or. vel_nu>1.07d0*d1) then
      vel_nu = d1
      call vel_matrices()
   end if
 end subroutine update_nu_r


!-------------------------------------------------------------------------
! add in terms from variable viscosity
!-------------------------------------------------------------------------
 subroutine update_nonlin_terms()
   use les
   use velocity
   implicit none
   type (coll) :: c1,c2,c3
   type (spec) :: s1,s2,s3
   type (phys) :: q(3), g(3), nu_
   logical, save :: set = .false.
   integer :: n, i, j
   _loop_km_vars
						   ! grad(nu_t)
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

   ! N'  = (1/Re).nu'.Lap(u) + nu_t.(-4\hat{z}) + grad(nu_t).2S
   ! nu' = Re*nu_t + 1 - nu_r
   nu_%Re = d_Re*nu%Re + (1d0-vel_nu)
   do i = 1, 3
      q(i)%Re = (1d0/d_Re) * nu_%Re * (-q(i)%Re)
      if(i==3) q(i)%Re = q(i)%Re + nu%Re*(-4d0) 
      do j = 1, 3
         q(i)%Re = q(i)%Re  +  g(j)%Re * S(j,i)%Re * 2d0
      end do
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
! write reference viscocity, max(nu_t), Re_tau
!-------------------------------------------------------------------------
 subroutine write_Retau()
   use les
   use io
   implicit none
   logical, save :: set=.false.

   if(mpi_rnk/=0) return

   if(.not.set) then
      set = .true.
      open(54,status='unknown',file='vel_Retau.dat')
   end if
   
   write(54,'(4e16.8)') tim_t, Re_tau, vel_nu, nu_max
   
 end subroutine write_Retau


!-------------------------------------------------------------------------
! (max(nu_t))(r)
!-------------------------------------------------------------------------
 subroutine save_nuprof()
   use les
   use io
   implicit none
   double precision :: d_(i_N), numax(i_N), numean(i_N)
   character(4) :: cnum
   integer :: n,n_

   numax = 0d0
   do n = 1, mes_D%pN
      n_ = mes_D%pNi + n - 1
      numax(n_) = maxval(nu%Re(:,:,n))
   end do
   numean = nuc%Re(:,0)

#ifdef _MPI
   call mpi_allreduce(numax, d_, i_N, mpi_double_precision,  &
      mpi_max, mpi_comm_world, mpi_er)
   numax = d_
#endif

   if(mpi_rnk/=0) return
   write(cnum,'(I4.4)') io_save1
   open(11, status='unknown', file='nu_prof'//cnum//'.dat')  
   do n = 1, i_N
      write(11,'(4e16.8)')  mes_D%r(n,1), numax(n), numean(n)
   end do
   close(11)
   
 end subroutine save_nuprof
