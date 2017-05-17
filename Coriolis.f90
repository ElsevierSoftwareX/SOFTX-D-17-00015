!*************************************************************************
! To activate the Coriolis term, the file  coriolis.in  must be in the 
! working directory.  It reads two parameters, the Ekman number Ek0,
! and the angle between the axis of the pipe and the rotation axis, phi.
!*************************************************************************
#define var_null util
#include "../program/main.f90"
!*************************************************************************

!-------------------------------------------------------------------------
! add coriolis force to nonlinear terms
!-------------------------------------------------------------------------
 subroutine util(FLAG)
   use io
   implicit none
   integer, intent(in) :: FLAG
   double precision, save :: Ek0, phi
   double precision, save :: Ro, Or(0:i_Th-1), Ot(0:i_Th-1), Oz, U(i_pN)
   logical, save :: fexist, set = .false.
   double precision :: d1
   type (phys) :: Nr,Nt,Nz
   type (spec) :: s1,s2,s3
   type (coll) :: c1,c2,c3
   integer :: n, m

   if(FLAG==1) then			! called just before write
      if(modulo(tim_step,i_save_rate1)==0) then
         call save_profiles()
      end if
   end if

   if(.not.set) then
      set = .true.
      inquire(file='coriolis.in', exist=fexist)
      if(.not.fexist) return
      if(mpi_rnk==0) print*, 'Including Coriolis force!'
      open(99,file='coriolis.in', status='old')
      read(99,*) Ek0, phi
      close(99)
      if(i_Mp/=1) stop 'add_coriolis: i_Mp/=1'
      d1 = 2d0*d_PI/dble(i_Th)
      do m = 0, i_Th-1
         Or(m) =  dsin(phi)*dcos(d1*m)	! \hat{Omeg}, unit vector
         Ot(m) = -dsin(phi)*dsin(d1*m)	! along axis of rotation
         Oz    =  dcos(phi)
      end do
      do n = 1, mes_D%pN
         U(n) = vel_U(mes_D%pNi-1+n)
      end do
   end if

   if(.not.fexist) return   

   if(FLAG/=2) return		! called after eval of nonlin terms

   Ro = d_Re * Ek0		! Ro = Re_0 Ek_0
   d1 = -1d0/(4d0*Ro)		! N += -1/(4Ro) (\hat{Omeg} x (U+u'))
   do n = 1, mes_D%pN
      do m = 0, i_Th-1
         Nr%Re(:,m,n) = d1*(Ot(m)*(vel_z%Re(:,m,n)+U(n)) - Oz*vel_t%Re(:,m,n))
         Nt%Re(:,m,n) = d1*(Oz*vel_r%Re(:,m,n) - Or(m)*(vel_z%Re(:,m,n)+U(n)))
         Nz%Re(:,m,n) = d1*(Or(m)*vel_t%Re(:,m,n) - Ot(m)*vel_r%Re(:,m,n))
      end do
   end do
   call tra_phys2spec(Nr, s1)
   call tra_phys2spec(Nt, s2)
   call tra_phys2spec(Nz, s3)
   call var_spec2coll(s1,c1, s2,c2, s3,c3)
   call var_coll_add(c1, vel_Nr)
   call var_coll_add(c2, vel_Nt)
   call var_coll_add(c3, vel_Nz)
   
 end subroutine util


!-------------------------------------------------------------------------
! write to files for non-newtonian data
!-------------------------------------------------------------------------
 subroutine save_profiles()
   use io
   implicit none
   double precision :: r(i_N), vz(-i_N:i_N), vz_(-i_N:i_N)
   character(4) :: cnum
   integer :: i,n, m1,m2,m3, pN,N1,N2

   r = mes_D%r(:,1)
   write(cnum,'(I4.4)') io_save1

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
