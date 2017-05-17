!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
 module velocity
!*************************************************************************
   use variables
   use transform
   use timestep
   implicit none
   save

   type (phys) :: vel_r
   type (phys) :: vel_t
   type (phys) :: vel_z
   type (phys) :: vel_curlr
   type (phys) :: vel_curlt
   type (phys) :: vel_curlz
   type (coll) :: vel_ur
   type (coll) :: vel_ut
   type (coll) :: vel_uz
   type (coll) :: vel_Nr
   type (coll) :: vel_Nt
   type (coll) :: vel_Nz
   double precision :: vel_nu
   double precision :: vel_Pr0
   double precision :: vel_U(i_N)
   double precision :: vel_Up(i_N)

   type (lumesh), private :: LDp(0:i_pH1), LDm(0:i_pH1)
   type (lumesh), private :: LDz(0:i_pH1), LNp(0:i_pH1)
   type (mesh),   private :: Ltp(0:i_pH1), Ltm(0:i_pH1)
   type (mesh),   private :: Ltz(0:i_pH1)
   type (coll),   private :: Nr_,Nt_,Nz_,ur_,ut_,uz_

   type (coll), private :: c1,c2,c3
   type (spec), private :: s1,s2,s3
   type (phys), private :: p

   
 contains

!------------------------------------------------------------------------
!  initialise velocity/pressure field
!------------------------------------------------------------------------
   subroutine vel_precompute()
      call var_coll_init(vel_ur)
      call var_coll_init(vel_ut)
      call var_coll_init(vel_uz)
      vel_nu  = 1d0
      vel_Pr0 = 0d0
      vel_U   = 1d0 - mes_D%r(:,2)
      vel_Up  = - 2d0*mes_D%r(:,1)
   end subroutine vel_precompute


!------------------------------------------------------------------------
!   symmetric, A(th)=A(-th);  s&r, A(th,z)=A(-th,z+L/2)
!------------------------------------------------------------------------
   subroutine vel_imposesym()   
      double precision :: a
      integer :: n
      _loop_km_vars

      call var_imposesym(0, vel_ur,vel_ut,vel_uz)
      if(b_mirrorsym)  &
         call var_imposesym(1, vel_ur,vel_ut,vel_uz)
      if(b_shiftrefl)  &
         call var_imposesym(2, vel_ur,vel_ut,vel_uz)
      if(b_shiftrott)  &
         call var_imposesym(3, vel_ur,vel_ut,vel_uz)

   end subroutine vel_imposesym 


!------------------------------------------------------------------------
!  convert 
!------------------------------------------------------------------------
   subroutine vel_rt2pm(r,t, up,um)
      type (coll), intent(in)  :: r,t
      type (coll), intent(out) :: up,um
      double precision :: rRe(i_N), rIm(i_N), tRe(i_N), tIm(i_N)
      _loop_km_vars
      _loop_km_begin
         rRe = r%Re(:,nh)
         rIm = r%Im(:,nh)
         tRe = t%Re(:,nh)
         tIm = t%Im(:,nh)
         up%Re(:,nh) = rRe - tIm
         up%Im(:,nh) = rIm + tRe
         um%Re(:,nh) = rRe + tIm
         um%Im(:,nh) = rIm - tRe
      _loop_km_end
   end subroutine vel_rt2pm


   subroutine vel_pm2rt(up,um, r,t)
      implicit none
      type (coll), intent(in)  :: up,um
      type (coll), intent(out) :: r,t
      double precision :: pRe(i_N),mRe(i_N),pIm(i_N),mIm(i_N)
      _loop_km_vars
      _loop_km_begin
         pRe = up%Re(:,nh)
         pIm = up%Im(:,nh)
         mRe = um%Re(:,nh)
         mIm = um%Im(:,nh)
         r%Re(:,nh) =  0.5d0 * (pRe + mRe)
         r%Im(:,nh) =  0.5d0 * (pIm + mIm)
         t%Re(:,nh) =  0.5d0 * (pIm - mIm)
         t%Im(:,nh) = -0.5d0 * (pRe - mRe)
      _loop_km_end
   end subroutine vel_pm2rt
   

!------------------------------------------------------------------------
!  precomputation of matrices for timestepping
!------------------------------------------------------------------------
   subroutine vel_matrices()
      double precision :: d1, d2
      logical, save :: set = .false.

			! pressure matrices
      if(.not.set) call tim_lumesh_init( 0,1,0d0,1d0, LNp)
      set = .true.

			! lhs matrices for uz and u+,u-
      d1 =  1d0/tim_dt
      d2 = -d_implicit*vel_nu/d_Re 
      call tim_lumesh_init( 1,0,d1,d2, LDp)
      call tim_lumesh_init(-1,0,d1,d2, LDm)
      call tim_lumesh_init( 0,0,d1,d2, LDz)

      			! timestepping matrices for rhs
      d1 =  1d0/tim_dt
      d2 =  (1d0-d_implicit)*vel_nu/d_Re
      call tim_mesh_init( 1,d1,d2, Ltp)
      call tim_mesh_init(-1,d1,d2, Ltm)
      call tim_mesh_init( 0,d1,d2, Ltz)

         		! get influence matrices
      call vel_adjustFlux(0)
      call vel_adjPPE(0)
      
   end subroutine vel_matrices
   

!------------------------------------------------------------------------
!  adjust flux if fixed; fix to 1/2, i.e. disturbance has zero mean flux
!------------------------------------------------------------------------
   subroutine vel_adjustFlux(F)
      integer, intent(in) :: F
      double precision, save :: Ui(i_N), d1,d2,d3
      integer :: info
      
      if(.not.b_const_flux) return
      if(mpi_rnk/=0) return

      if(F==0) then
         Ui(:)   = 1d0
         Ui(1)   = 0d0
         Ui(i_N) = 0d0
         call dgbtrs('N', i_N, i_KL, i_KL, 1, LDz(0)%M, 3*i_KL+1,  &
                     LDz(0)%ipiv, Ui, i_N, info)
         if(info/=0) stop 'vel_adjustFlux: err 1'
         d1 = 2d0*dot_product(Ui, mes_D%intrdr)
         if(d1==0d0) stop 'vel_adjustFlux: err 2'
      
      else if(F==1) then
         d2 = 2d0*dot_product(vel_uz%Re(:,0), mes_D%intrdr)
         d3 = -d2/d1
         vel_uz%Re(:,0) = vel_uz%Re(:,0) + d3*Ui
         vel_Pr0 = vel_Pr0 + d3*d_Re/4d0
      end if

   end subroutine vel_adjustFlux


!------------------------------------------------------------------------
!  PPE formulation + influence matrix adjusts to correct bcs.
!------------------------------------------------------------------------
   subroutine vel_adjPPE(F)
      integer, intent(in) :: F
      double precision, save :: A(4,4,0:i_pH1)
      double precision, allocatable, save :: U(:,:,:) !i_N,0:i_pH1,6
      double precision :: BRe(4,0:i_pH1),BIm(4,0:i_pH1), aR(4),aI(4)
      integer :: j
      _loop_km_vars
      
      if(F==0) then                             ! precompute, get U
         if(.not.allocated(U))  &
            allocate( U(i_N,0:i_pH1,6) )
         do j = 1, 4
            call var_coll_init(c1)
            call var_coll_init(c2)
            call var_coll_init(c3)
            if(j==1) then
               c1%Re(i_N,:) = 1d0
               call tim_lumesh_invert(0,LDp, c1)
               U(:,:,1) = c1%Re
            else if(j==2) then
               c2%Re(i_N,:) = 1d0
               call tim_lumesh_invert(0,LDm, c2)
               U(:,:,2) = c2%Re
            else if(j==3) then
               c3%Im(i_N,:) = 1d0
               call tim_lumesh_invert(0,LDz, c3)
               U(:,:,3) = c3%Im
            else
               c1%Re(i_N,:) = -1d0
               call tim_lumesh_invert(1,LNp, c1)
               call var_coll_grad(c1, c1,c2,c3)
               call vel_rt2pm(c1,c2, c1,c2)
               U(:,:,4) = c1%Re
               U(:,:,5) = c2%Re
               U(:,:,6) = c3%Im
            end if
            call vel_evalBC(c1,c2,c3, BRe,BIm)
            A(:,j,:) =  BRe
         end do
         _loop_km_begin
            if(k==0 .and. m==0) cycle
            call mes_mat_invert(4,A(1,1,nh),4)
         _loop_km_end

      else if(F==1) then			! get p, project RHS
         call vel_pm2rt(vel_ur,vel_ut, c1,c2)
         call var_coll_div(c1,c2,vel_uz, c1)
         call tim_zerobc(c1)
         call tim_lumesh_invert(1,LNp, c1)
         call var_coll_grad(c1, c1,c2,c3)
         call vel_rt2pm(c1,c2, c1,c2)
         call var_coll_sub(c1, vel_ur)
         call var_coll_sub(c2, vel_ut)
         call var_coll_sub(c3, vel_uz)

      else if(F==2) then                        ! correct the BCs
         call vel_evalBC(vel_ur,vel_ut,vel_uz, BRe,BIm)
         _loop_km_begin
            if(k==0 .and. m==0) cycle
            aR = -matmul(A(:,:,nh),BRe(:,nh))
            aI = -matmul(A(:,:,nh),BIm(:,nh))
            vel_ur%Re(:,nh) = vel_ur%Re(:,nh) + aR(1)*U(:,nh,1)
            vel_ur%Im(:,nh) = vel_ur%Im(:,nh) + aI(1)*U(:,nh,1)
            vel_ut%Re(:,nh) = vel_ut%Re(:,nh) + aR(2)*U(:,nh,2)
            vel_ut%Im(:,nh) = vel_ut%Im(:,nh) + aI(2)*U(:,nh,2)
            vel_uz%Re(:,nh) = vel_uz%Re(:,nh) - aI(3)*U(:,nh,3)
            vel_uz%Im(:,nh) = vel_uz%Im(:,nh) + aR(3)*U(:,nh,3)
            vel_ur%Re(:,nh) = vel_ur%Re(:,nh) + aR(4)*U(:,nh,4)
            vel_ur%Im(:,nh) = vel_ur%Im(:,nh) + aI(4)*U(:,nh,4)
            vel_ut%Re(:,nh) = vel_ut%Re(:,nh) + aR(4)*U(:,nh,5)
            vel_ut%Im(:,nh) = vel_ut%Im(:,nh) + aI(4)*U(:,nh,5)
            vel_uz%Re(:,nh) = vel_uz%Re(:,nh) - aI(4)*U(:,nh,6)
            vel_uz%Im(:,nh) = vel_uz%Im(:,nh) + aR(4)*U(:,nh,6)
         _loop_km_end
      end if

   end subroutine vel_adjPPE

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   subroutine vel_evalBC(up,um,uz, BRe,BIm)
      type (coll),      intent(in)  :: up,um,uz
      double precision, intent(out) :: BRe(4,0:i_pH1), BIm(4,0:i_pH1)
      double precision :: drRe,drIm, urRe,urIm, utRe,utIm, uzRe,uzIm
      double precision :: d, s, a_, b(0:i_M*i_Mp)
      integer :: n,m_
      _loop_km_vars
      
      _loop_km_begin
         if(k==0 .and. m==0) cycle
         m_ = m*i_Mp
         a_ = d_alpha*k
         n = i_N

         urRe =  0.5d0 * (up%Re(n,nh) + um%Re(n,nh))
         urIm =  0.5d0 * (up%Im(n,nh) + um%Im(n,nh))
         utRe =  0.5d0 * (up%Im(n,nh) - um%Im(n,nh))
         utIm = -0.5d0 * (up%Re(n,nh) - um%Re(n,nh))
         uzRe =  uz%Re(n,nh)
         uzIm =  uz%Im(n,nh)
         drRe =  0.5d0 * dot_product(  &
            mes_D%dr1(:,1), up%Re(i_N-i_KL:i_N,nh)+um%Re(i_N-i_KL:i_N,nh))
         drIm =  0.5d0 * dot_product(  &
            mes_D%dr1(:,1), up%Im(i_N-i_KL:i_N,nh)+um%Im(i_N-i_KL:i_N,nh))
            
         BRe(1,nh) = up%Re(n,nh)
         BIm(1,nh) = up%Im(n,nh)
         BRe(2,nh) = um%Re(n,nh)
         BIm(2,nh) = um%Im(n,nh)
         BRe(3,nh) = -uzIm
         BIm(3,nh) =  uzRe
         d = mes_D%r(n,-1)
         BRe(4,nh) = urRe*d + drRe - d*m_*utIm - a_*uzIm
         BIm(4,nh) = urIm*d + drIm + d*m_*utRe + a_*uzRe
      _loop_km_end

   end subroutine vel_evalBC


!------------------------------------------------------------------------
!  Evaluate in physical space  u  and  curl(u)
!------------------------------------------------------------------------
   subroutine vel_transform()
      
      call var_coll2spec(vel_ur,s1, c2=vel_ut,s2=s2, c3=vel_uz,s3=s3)
      call tra_spec2phys( s1, vel_r)
      call tra_spec2phys( s2, vel_t)
      call tra_spec2phys( s3, vel_z)

      call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
      call var_coll2spec(c1,s1, c2=c2,s2=s2, c3=c3,s3=s3)
      call tra_spec2phys( s1, vel_curlr)
      call tra_spec2phys( s2, vel_curlt)
      call tra_spec2phys( s3, vel_curlz)

   end subroutine vel_transform


!-------------------------------------------------------------------------
!  nonlinear terms for velocity
!------------------------------------------------------------------------
   subroutine vel_nonlinear()
         			! advection  u x curlu
      p%Re = vel_t%Re*vel_curlz%Re - vel_z%Re*vel_curlt%Re
      call tra_phys2spec(p, s1)
      p%Re = vel_z%Re*vel_curlr%Re - vel_r%Re*vel_curlz%Re
      call tra_phys2spec(p, s2)
      p%Re = vel_r%Re*vel_curlt%Re - vel_t%Re*vel_curlr%Re
      call tra_phys2spec(p, s3)
      call var_spec2coll(s1,vel_Nr, s2=s2,c2=vel_Nt, s3=s3,c3=vel_Nz)
      
      call vel_addHPF()

   end subroutine vel_nonlinear

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine vel_addHPF()
      double precision :: a(i_N), b(i_N)
      _loop_km_vars
                  			! force from background HPF
      b = -vel_Up
      _loop_km_begin
         a = d_alpha*k * vel_U
         vel_Nr%Re(:,nh) = vel_Nr%Re(:,nh) + a*vel_ur%Im(:,nh)
         vel_Nr%Im(:,nh) = vel_Nr%Im(:,nh) - a*vel_ur%Re(:,nh)
         vel_Nt%Re(:,nh) = vel_Nt%Re(:,nh) + a*vel_ut%Im(:,nh)
         vel_Nt%Im(:,nh) = vel_Nt%Im(:,nh) - a*vel_ut%Re(:,nh)
         vel_Nz%Re(:,nh) = vel_Nz%Re(:,nh) + a*vel_uz%Im(:,nh)  &
                                           + b*vel_ur%Re(:,nh) 
         vel_Nz%Im(:,nh) = vel_Nz%Im(:,nh) - a*vel_uz%Re(:,nh)  & 
                                           + b*vel_ur%Im(:,nh)
      _loop_km_end
       					! additional pressure if fixed flx
      if(b_const_flux .and. mpi_rnk==0) then
         vel_Pr0 = dot_product(vel_uz%Re(i_N-i_KL:,0),mes_D%dr1(:,1))
         vel_Pr0 = vel_Pr0/(-2d0)
         vel_Nz%Re(:,0) = vel_Nz%Re(:,0) + 4d0*vel_Pr0/d_Re
      end if

   end subroutine vel_addHPF
   

!------------------------------------------------------------------------
!  Advance equations one timestep.
!------------------------------------------------------------------------
   subroutine vel_step()
				     	! get rhs = A u_ + N
      call vel_rt2pm(vel_Nr,vel_Nt, c1,c2)
      call vel_rt2pm(ur_,ut_, vel_ur,vel_ut)
      call tim_meshmult(1,Ltp,vel_ur,c1,  vel_ur)
      call tim_meshmult(1,Ltm,vel_ut,c2,  vel_ut)
      call tim_meshmult(0,Ltz,uz_,vel_Nz, vel_uz)
      call vel_adjPPE(1)
        				! invert
      call tim_zerobc(vel_ur)
      call tim_zerobc(vel_ut)
      call tim_zerobc(vel_uz)
      call tim_lumesh_invert(0,LDp, vel_ur)
      call tim_lumesh_invert(0,LDm, vel_ut)
      call tim_lumesh_invert(0,LDz, vel_uz)
      call vel_adjPPE(2)
      call vel_pm2rt(vel_ur,vel_ut, vel_ur,vel_ut)
         				! 
      call vel_adjustFlux(1)

      if(mpi_rnk==0) then
         vel_ur%Re(:,0) = 0d0
         vel_ur%Im(:,0) = 0d0
         vel_ut%Im(:,0) = 0d0
         vel_uz%Im(:,0) = 0d0
      end if
            
   end subroutine vel_step
   
   
!------------------------------------------------------------------------
!  predictor with euler nonlinear terms
!------------------------------------------------------------------------
   subroutine vel_predictor()

      call var_coll_copy(vel_ur, ur_)
      call var_coll_copy(vel_ut, ut_)
      call var_coll_copy(vel_uz, uz_)
      call var_coll_copy(vel_Nr, Nr_)
      call var_coll_copy(vel_Nt, Nt_)
      call var_coll_copy(vel_Nz, Nz_)
      call vel_step()
      
   end subroutine vel_predictor


!------------------------------------------------------------------------
!  corrector iteration with Crank-Nicolson non-lin term
!------------------------------------------------------------------------
   subroutine vel_corrector()
      type (coll) :: r,t,z

      call var_coll_copy(vel_ur, r)
      call var_coll_copy(vel_ut, t)
      call var_coll_copy(vel_uz, z)
      call tim_nlincorr(Nr_, vel_Nr)
      call tim_nlincorr(Nt_, vel_Nt)
      call tim_nlincorr(Nz_, vel_Nz)
      call vel_step()

      call var_coll_sub(vel_ur, r)
      call var_coll_sub(vel_ut, t)
      call var_coll_sub(vel_uz, z)
      call tim_measurecorr(r,t,z)

   end subroutine vel_corrector

   
!------------------------------------------------------------------------
!  get cfl max dt due to flow field
!------------------------------------------------------------------------
   subroutine vel_maxtstep()
      double precision :: d,mx, dt(3),dt_(3), r_(0:i_N+1)
      integer :: N_,n, n__

      N_ = mes_D%pN
      do n = 0, N_+1
         n__ = n+mes_D%pNi-1
         if(n__>=1 .and. n__<=mes_D%N)  r_(n) = mes_D%r(n__,1)
      end do

      dt = 1d99

      do n = 1, N_

         if(n==1 .and. mpi_rnk==0) then
            d = r_(2) - r_(1)
         else if(n==N_ .and. mpi_rnk==mpi_sze-1) then
            d = r_(N_) - r_(N_-1)
         else
            d = min( r_(n)-r_(n-1), r_(n+1)-r_(n) )
         end if
         mx = maxval( dabs(vel_r%Re(:,:,n)) )
         if(mx/=0d0) dt(1) = min( dt(1), d/mx )

         d = 2d0*d_PI/dble(i_Th*i_Mp) 		!---  *r_(n)? ---
         mx = maxval( dabs(vel_t%Re(:,:,n)) )
         if(mx/=0d0) dt(2) = min( dt(2), d/mx )

         d = 2d0*d_PI/(d_alpha*i_Z)
         mx = maxval( dabs(vel_z%Re(:,:,n) + (1d0-r_(n)*r_(n))))
         if(mx/=0d0) dt(3) = min( dt(3), d/mx )
      end do

#ifdef _MPI
      call mpi_allreduce(dt, dt_, 3, mpi_double_precision,  &
         mpi_min, mpi_comm_world, mpi_er)
      dt = dt_
#endif
      tim_cfl_dt = minval(dt)
      if(tim_cfl_dt==dt(1)) tim_cfl_dir=1
      if(tim_cfl_dt==dt(2)) tim_cfl_dir=2
      if(tim_cfl_dt==dt(3)) tim_cfl_dir=3
      
   end subroutine vel_maxtstep


!*************************************************************************
 end module velocity
!*************************************************************************
