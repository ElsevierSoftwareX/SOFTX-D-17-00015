!*************************************************************************
!  Newton vector:
!    x(1) = theta shift
!    x(2) = z shift
!    x(3) = period of orbit
!    x(4) = continuation parameter
!    x(9:)= velocity vector
!  Constraints:
!    x . dx/dth = 0,  x . dx/dz
!       no update parallel to theta or z shift.
!    (F(x)-x). dx/dt = 0 .
!       no update along direction of trajectory.
!  Jacobian approximation:
!    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
!
!  File shifts.in :
!     Ndt	  	number of timesteps taken in one period
!     Sth Sz T		initial guesses for shift and period
!     
!*************************************************************************
 module orbit
   use parameters
   implicit none
   save

   logical          :: T_fixed
   integer          :: m, nits, ms, ncgd
   double precision :: rel_err
   double precision :: del, mndl, mxdl
   double precision :: gtol, epsJ

   logical          :: Ctnu
   integer          :: ctntmn, ctntmx, ctp
   double precision :: ctdp0

   integer, parameter :: n = 8 + 6*i_N*(i_pH1+1)
   double precision   :: tol, scaleU(3), scaleT, ctpsc, ctstp(n), ctx_(n)
   integer :: info, ndts(0:9)
   real    :: d_start, d_stop

 end module orbit


!*************************************************************************
#include "../parallel.h"
#include "NewtonHook.f90"
#include "GMRESm.f90"
 PROGRAM MAIN
!*************************************************************************
   use velocity
   use io
   use orbit
   use newton
   implicit none
   double precision :: d
   external :: getrhs, multJ, multJp, saveorbit 
   double precision, external :: dotprod, dotprod_ms
   integer :: nms, i1, i2, i3, p, cti
   character(4) :: cnum

   open(99,status='old',file='params.in')
   read(99,*) i1
   read(99,*) m, nits, ms, ncgd
   read(99,*) rel_err
   read(99,*) del, mndl, mxdl
   read(99,*) gtol, epsJ
   close(99)
   T_fixed = (i1==1)
   nms = n * ms
   allocate(new_x(nms))
   allocate(new_fx(nms))
  
   inquire(file='params_ct.in', exist=Ctnu)
   if(Ctnu) then
      open(99,status='old',file='params_ct.in')
      read(99,*) ctntmn, ctntmx, ctp
      read(99,*) ctdp0
   end if

   open(99,status='old',file='shifts.in')
   do p = 0, ms-1
      read(99,*) ndts(p)
   end do
   do p = 0, ms-1
      new_x(p*n+1:p*n+8) = 0d0
      read(99,*) new_x(p*n+1), new_x(p*n+2), new_x(p*n+3)
   end do
   close(99)

   call initialise()  
   call vel_imposesym()
   call get_scales()
   call rtz2vec(n,vel_ur,vel_ut,vel_uz, new_x)
   tim_t = 0d0
   call load_ms_states(0)
   call search_sz()

   d = dotprod_ms(-1,new_x,new_x)
   tol  = rel_err * dsqrt(d)
   del  = del  * dsqrt(d)
   mndl = mndl * dsqrt(d)
   mxdl = mxdl * dsqrt(d)

   if(Ctnu) then
      ctx_ = new_x
      call load_ms_states(1)
      ctstp = new_x - ctx_
      if(maxval(dabs(ctstp))==0d0)  &
         ctstp(4::n) = ctstp(4::n) + 1d0
      ctx_  = new_x
      new_x = new_x + ctstp*ctdp0
      if(ctp==1) ctpsc = d_Re    / dsqrt(d)
      if(ctp==2) ctpsc = d_alpha / dsqrt(d)
      new_x(4::n) = new_x(4::n)/ctpsc
      ctstp(4::n) = ctstp(4::n)/ctpsc
      ctx_ (4::n) = ctx_ (4::n)/ctpsc
   end if

   cti = 0
   do while(.true.) 
      info = 1
      if(mpi_rnk/=0) info = 0
      call newtonhook(getrhs, multJ, multJp, saveorbit, dotprod_ms, &
                      m, nms, gtol, tol, del, mndl, mxdl, nits, info)
      if(ncgd>0 .and. info==0) call getEigen() 
      call io_closefiles()
      if(.not.Ctnu) exit
      cti = cti + 1
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      if(mpi_rnk==0) then
         write(cnum,'(I4)') cti
         if(cti<1000) cnum(1:1) = '0'
         if(cti< 100) cnum(2:2) = '0'
         if(cti<  10) cnum(3:3) = '0'
         call system('sleep 20s ; mkdir '//cnum//' ; sleep 20s ;'// &
            'mv *.dat '//cnum//' ; sleep 20s')
         d = sum(new_x(2::n))/sum(new_x(3::n))  !c
         open(99,status='unknown',access='append',file='newton.dat_ct')   
         write(99,'(2I5,3e20.12)') cti, new_nits, d_Re, d_alpha, d
         close(99)   
      end if
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      if(info/=0) exit 
      call var_precompute()
      call tim_precompute()
      call vel_precompute()
      call  io_precompute()
      call  io_openfiles()
      d = 1.1d0
      if(new_nits>ctntmx) d=0.5d0
      if(new_nits<ctntmn) d=2.0d0
      ctstp = new_x - ctx_
      ctx_  = new_x
      new_x = new_x + ctstp*d
      if(T_fixed) then
         d = 1d0
         if(new_gits>(2*m)/3) d=2d0
         if(new_gits<10) d=0.5d0
         ndts = int(dble(ndts)*d+0.5d0)
         do p = 1, 3
            new_x(p::n) = new_x(p::n)*d
            ctx_ (p::n) = ctx_ (p::n)*d
         end do
      end if
   end do

   call cleanup()      
   stop

 contains

!-------------------------------------------------------------------------
!  Initialisation
!-------------------------------------------------------------------------
   subroutine initialise()
      logical :: file_exist

      call mpi_precompute()
      if(mpi_rnk==0) then
         call system('touch PRECOMPUTING')
         call system('echo $HOSTNAME > HOST')
      end if

      if(mpi_rnk==0)  print*, 'precomputing function requisites...'
      call par_precompute()
      call mes_precompute()
      call var_precompute()
      call tra_precompute()
      call tim_precompute()
      call vel_precompute()
      call  io_precompute()
   
      if(mpi_rnk==0)  print*, 'loading state...'
      tim_dt = 1d99
      call io_load_state()
      call vel_matrices()

      if(mpi_rnk==0)  print*, 'initialising output files...'
      call io_openfiles()

      if(mpi_rnk==0) then
         open (99, file='PRECOMPUTING')
         close(99, status='delete')
         open (99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
         close(99)
         print*, 'timestepping.....'
      end if
      
      call clk_time(d_start)
   
   end subroutine initialise

!-------------------------------------------------------------------------
!  Program closure
!-------------------------------------------------------------------------
   subroutine cleanup()
      logical :: file_exist
   
      if(mpi_rnk==0) then
         print*, 'cleanup...'
         call clk_time(d_stop)
         print*, ' sec/step  = ', (d_stop-d_start)/real(tim_step)
#ifdef _CPUTIME
         print*, ' CPU time  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
         print*, ' WALL time = ', int((d_stop-d_start)/6d1), ' mins.'
#endif 
      end if
      
!      call io_save_state()
!      call io_closefiles()

      if(mpi_rnk==0) then
         inquire(file='RUNNING', exist=file_exist)
         if(file_exist) open(99, file='RUNNING')
         if(file_exist) close(99, status='delete')
      end if      

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif
      if(mpi_rnk==0) print*, '...done!'

   end subroutine cleanup


!-------------------------------------------------------------------------
   subroutine clk_time(t)
      real, intent(out) :: t
#ifdef _CPUTIME
      call cpu_time(t)
#else
      integer, save :: ct,ctrt,ctmx,ct_=0,wrap=0
      call system_clock(ct,ctrt,ctmx)
      if(ct<ct_) wrap = wrap + 1
      ct_ = ct
      t = (real(ct)+real(ctmx)*wrap)/real(ctrt)
#endif
   end subroutine clk_time
         
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

 subroutine get_scales()
   use velocity
   use newton
   use orbit
   implicit none
   double precision :: d,Er,Et,Ez,Ez00, Ek(0:i_K1), Em(0:i_M1)
   double precision, external :: dotprod
   type (coll) :: c1
   integer :: p
   
   call var_coll_norm(vel_ur, Er,Ek,Em)
   call var_coll_norm(vel_ut, Et,Ek,Em)
   call var_coll_norm(vel_uz, Ez,Ek,Em)
   call var_coll_copy(vel_uz, c1)
   if(mpi_rnk==0) c1%Re(:,0) = 0d0
   call var_coll_norm(c1, Ez00,Ek,Em)
   Ez00 = Ez - Ez00
   Ez = Ez - Ez00

   scaleU(1) = dsqrt(0.2d0*Ez00/Er)
   scaleU(2) = dsqrt(0.2d0*Ez00/Et)
   scaleU(3) = dsqrt(0.5d0*Ez00/Ez)

   call rtz2vec(n,vel_ur,vel_ut,vel_uz, new_x)
   d = dotprod(-1,new_x,new_x)
   scaleT   = new_x(3) / dsqrt(d)
   do p = 0, ms-1
      new_x(p*n+1:p*n+3) = new_x(p*n+1:p*n+3) / scaleT
   end do

   if(mpi_rnk==0)  &
      print*, 'Scales:', real(scaleU), real(scaleT)
   
 end subroutine get_scales
 

!-------------------------------------------------------------------------
! get shooting points
!-------------------------------------------------------------------------
 subroutine load_ms_states(i)
   use newton
   use orbit
   use io
   implicit none
   integer, intent(in) :: i
   double precision :: t_, x(n), y(n)
   character(2) :: cnum
   integer :: p, p1, e,f
   logical :: fexist

   if(Ctnu) then
      if(i==0) e=nf90_open('state.cdf.in',nf90_nowrite, f)
      if(i==1) e=nf90_open('state10.cdf.in',nf90_nowrite, f)
      if(e/=nf90_noerr) return
      if(ctp==1) e=nf90_get_att(f,nf90_global,'Re', t_)
      if(ctp==2) e=nf90_get_att(f,nf90_global,'alpha', t_)
      e=nf90_close(f)
      new_x(4::n) = t_
   end if

   do p = 1-i, ms-1
      write(cnum, '(I2)') p + i*10
      io_statefile = 'state'//cnum(2-i:2)//'.cdf.in'
      inquire(file='state'//cnum(2-i:2)//'.cdf.in', exist=fexist)
      if(fexist) then 
         t_ = tim_t
         call io_load_state()
         call vel_imposesym()
         call rtz2vec(n,vel_ur,vel_ut,vel_uz, new_x(p*n+1))
         tim_t = t_
      else
         p1 = p-1
         x = new_x(p1*n+1:p1*n+n)
         call steporbit(ndts(p1),x, y)
         new_x(p*n+9:p*n+n) = y(9:)         
      end if
   end do
   
 end subroutine load_ms_states


!-------------------------------------------------------------------------
! if st=0 and sz=-1, search for initial sz, for last ms only.
!-------------------------------------------------------------------------
 subroutine search_sz()
   use newton
   use orbit
   use io
   implicit none
   double precision, external :: dotprod
   double precision :: d1,d2,sz,sz_,st, x(n),y(n)
   type(coll) :: r,t,z,r_,t_,z_
   integer :: p,p1
   do p1 = 0, ms-1
      p = modulo(p1+1,ms)
      st = new_x(p1*n+1)*scaleT
      sz = new_x(p1*n+2)*scaleT
      if(dabs(st)>1d-8) cycle
      if(dabs(sz+1d0)>1d-8) cycle
      if(mpi_rnk==0) print*, 'search_sz: p1 = ', p1
      if(Ctnu) STOP 'search_sz: initial sz required for continuation.'
      x = new_x(p*n+1:p*n+n)
      y = new_x(p1*n+1:p1*n+n)
      y(2) = 0d0
      call steporbit(ndts(p1),y, y)
      call vec2rtz(n,y, r,t,z)
      d1 = 1d99
      sz = 0d0
      do while(sz<2d0*d_PI/d_alpha)
         call var_coll_shift(0d0,-sz,r, r_)
         call var_coll_shift(0d0,-sz,t, t_)
         call var_coll_shift(0d0,-sz,z, z_)
         call rtz2vec(n,r_,t_,z_, y)
         y = y - x
         d2 = dotprod(-1,y,y)
         if(d2<d1) sz_ = sz
         if(d2<d1) d1  = d2
         sz = sz + 0.01d0
      end do
      new_x(p1*n+2) = sz_/scaleT
      if(mpi_rnk==0) print*, 'search_sz: sz = ', real(sz_)
   end do
 end subroutine search_sz

!-------------------------------------------------------------------------
!  function to be minimised   
!-------------------------------------------------------------------------
 subroutine getrhs(nms,x, y)
   use parameters
   use orbit
   implicit none
   integer,          intent(in)  :: nms
   double precision, intent(in)  :: x(nms)
   double precision, intent(out) :: y(nms)
   double precision :: x_(n), y_(n*ms)
   integer :: p, p1
   
   do p = 0, ms-1
      p1 = modulo(p+1,ms)
      x_ = x(p*n+1:p*n+n)
      call steporbit(ndts(p),x_, y_(p1*n+1))
   end do
   y = y_ - x					! diff
   do p = 0, ms-1
      y(p*n+1:p*n+8) = 0d0			! constraints, rhs=0
   end do

 end subroutine getrhs


!-------------------------------------------------------------------------
!  Jacobian of function + lhs of constraints on update
!-------------------------------------------------------------------------
 subroutine multJ(nms,x, y)
   use parameters
   use newton
   use orbit,    only : n, ms, T_fixed, epsJ, Ctnu, ctstp
   use timestep, only : tim_dt
   implicit none
   integer,          intent(in)  :: nms
   double precision, intent(in)  :: x(nms)
   double precision, intent(out) :: y(nms)   
   double precision, external :: dotprod, dotprod_ms
   double precision :: eps, s(n*ms)
   integer :: p

    				! (F(x0+eps.x)-F(x0))/eps
   eps = dsqrt(dotprod_ms(1,x,x))
   if(eps==0d0)  stop 'multJ: eps=0 (1)'
   eps = epsJ * dsqrt(dotprod_ms(1,new_x,new_x)) / eps
   if(eps==0d0)  stop 'multJ: eps=0 (2)'
   y = new_x + eps*x
   call getrhs(nms,y, s)
   y = (s - new_fx) / eps

   do p = 0, ms-1
       				! no update in shift directions 
      call getshiftdir(1,n,new_x(p*n+1), s) 
      y(p*n+1) = dotprod(-1,s,x(p*n+1))
      call getshiftdir(2,n,new_x(p*n+1), s) 
      y(p*n+2) = dotprod(-1,s,x(p*n+1))
      				! no update in trajectory direction
      call steporbit(1,new_x(p*n+1), s)
      s(:n) = (s(:n) - new_x(p*n+1:p*n+n)) / tim_dt
      y(p*n+3) = dotprod(-1,s,x(p*n+1))
      				! special cases
      if(b_mirrorsym .or.  &
         b_shiftrefl )  y(p*n+1) = x(p*n+1)
      if(T_fixed)       y(p*n+3) = x(p*n+3)
   end do
   
   if(Ctnu)  y(4::n) = dotprod_ms(1,ctstp,x)

 end subroutine multJ
 

!-------------------------------------------------------------------------
!  preconditioner for multJ   
!-------------------------------------------------------------------------
 subroutine multJp(n, x)
   implicit none
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
 end subroutine multJp


!-------------------------------------------------------------------------
!  called at each newton iteration   
!-------------------------------------------------------------------------
 subroutine saveorbit()
   use newton
   use orbit
   use io
   implicit none
   double precision, external :: dotprod_ms
   double precision :: norm_x, p(3), normE_x, normE_fx
   integer :: p_, ndt
     
   ndt = sum(ndts(0:ms-1))
   norm_x = dsqrt(dotprod_ms(-1,new_x,new_x))
   call normE(new_x, normE_x)
   call normE(new_fx, normE_fx)
   
   if(mpi_rnk==0) then
      open(99,status='unknown',access='append',file='newton.dat')
      if(new_nits==0)  write(99,*) ndt, m, 3+6*i_N*(i_H1+1)
      write(99,'(2I6,4e13.5)')  &
         new_nits, new_gits, new_tol, new_del, new_tol/norm_x, norm_x
      close(99)

      open(99,status='unknown',access='append',file='newtonE.dat')
      write(99,'(I6,4e13.5)')  &
         new_nits, normE_fx/normE_x, normE_x
      close(99)

      open(99,status='unknown',access='append',file='shifts.dat')
      do p_ = 0, ms-1
         if(new_nits==0)  write(99,*) ndts(p_)
      end do
      do p_ = 0, ms-1
         p = new_x(p_*n+1:p_*n+3)*scaleT
         write(99,'(1I6,3e24.16)')  new_nits, p(1), p(2), p(3)
      end do
      close(99)
   end if

   do p_ = 0, ms-1
      call vec2rtz(n,new_x(p_*n+1), vel_ur,vel_ut,vel_uz)
      if(ms==1) io_save1 = new_nits
      if(ms/=1) io_save1 = new_nits*10 + p_
      call io_save_state()
      if(p_==0) call io_save_spectrum()
      if(p_==0) call io_save_meanprof()   
   end do  
   
 end subroutine saveorbit
 
 
!-------------------------------------------------------------------------
!  direction in phase space of a shift of x in theta or z 
!-------------------------------------------------------------------------
 subroutine getshiftdir(i,n,x, s)
   use variables
   implicit none
   integer,          intent(in)  :: i, n
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: s(n)
   double precision :: ds, d(i_N)
   type (coll) :: r,t,z
   _loop_km_vars
   
   call vec2rtz(n,x, r,t,z)

   _loop_km_begin
      if(i==1)  ds = dble(i_Mp*m)
      if(i==2)  ds = d_alpha*k
      d(:)       = -r%Im(:,nh)*ds
      r%Im(:,nh) =  r%Re(:,nh)*ds 
      r%Re(:,nh) =  d
      d(:)       = -t%Im(:,nh)*ds
      t%Im(:,nh) =  t%Re(:,nh)*ds 
      t%Re(:,nh) =  d
      d(:)       = -z%Im(:,nh)*ds
      z%Im(:,nh) =  z%Re(:,nh)*ds 
      z%Re(:,nh) =  d
   _loop_km_end

   call rtz2vec(n,r,t,z, s)
 
 end subroutine getshiftdir


!-------------------------------------------------------------------------
!  save state as a single vector
!-------------------------------------------------------------------------
 subroutine rtz2vec(n,r,t,z, v)
   use variables
   use orbit, only : scaleU
   implicit none
   integer,          intent(in)  :: n
   type (coll),      intent(in)  :: r,t,z
   double precision, intent(out) :: v(n)
   double precision, save :: d(i_N), s(3), w(0:i_pH1)
   logical, save :: set = .false.
   integer :: i, n_
   _loop_km_vars

   if(.not.set) then
      set = .true.
      d = dsqrt(dabs(mes_D%intrdr))
      _loop_km_begin
         w(nh) = 4d0/(4d0+(d_alpha*abs(k))+dble(i_Mp*abs(m))) 
      _loop_km_end
   end if

   i = 9
   _loop_km_begin
     do n_ = 1, i_N
        s = d(n_)*w(nh)*scaleU
        if(k==0 .and. m==0) s(3) = d(n_)
        v(i  ) = r%Re(n_,nh) * s(1)
        v(i+1) = r%Im(n_,nh) * s(1)
        v(i+2) = t%Re(n_,nh) * s(2)
        v(i+3) = t%Im(n_,nh) * s(2)
        v(i+4) = z%Re(n_,nh) * s(3)
        v(i+5) = z%Im(n_,nh) * s(3)
        i = i + 6
     end do
   _loop_km_end   
   
 end subroutine rtz2vec


!-------------------------------------------------------------------------
!  get state from vector
!-------------------------------------------------------------------------
 subroutine vec2rtz(n,v, r,t,z)
   use variables
   use orbit, only : scaleU
   implicit none
   integer,          intent(in)  :: n
   double precision, intent(in)  :: v(n)
   type (coll),      intent(out) :: r,t,z
   double precision, save :: d(i_N), s(3), w(0:i_pH1)
   logical, save :: set = .false.
   integer :: i, n_
   _loop_km_vars

   if(.not.set) then
      set = .true.
      d = dsqrt(dabs(mes_D%intrdr))
      _loop_km_begin
         w(nh) = 4d0/(4d0+(d_alpha*abs(k))+dble(i_Mp*abs(m)))
      _loop_km_end
   end if

   i = 9
   _loop_km_begin
     do n_ = 1, i_N
        s = d(n_)*w(nh)*scaleU
        if(k==0 .and. m==0) s(3) = d(n_)
        r%Re(n_,nh) = v(i  ) / s(1)
        r%Im(n_,nh) = v(i+1) / s(1)
        t%Re(n_,nh) = v(i+2) / s(2)
        t%Im(n_,nh) = v(i+3) / s(2)
        z%Re(n_,nh) = v(i+4) / s(3)
        z%Im(n_,nh) = v(i+5) / s(3)
        i = i + 6
     end do
   _loop_km_end   
   
 end subroutine vec2rtz


!-------------------------------------------------------------------------
! dot prod of line vector
!-------------------------------------------------------------------------
 double precision function dotprod(n_,a,b)
   use variables
   use orbit
   implicit none
   integer,          intent(in) :: n_
   double precision, intent(in) :: a(n), b(n)
   double precision :: d,d_
   integer :: n1,nn
   nn = 8 + 6*i_N*(var_H%pH1+1)
   n1 = 1
   if(mpi_rnk/=0 .or. n_==-1) n1 = 9
   d = dot_product(a(n1:nn),b(n1:nn))
#ifdef _MPI
   call mpi_allreduce(d, d_, 1, mpi_double_precision,  &
      mpi_sum, mpi_comm_world, mpi_er)
   d = d_
#endif
   dotprod = d
 end function dotprod


!-------------------------------------------------------------------------
! dot prod for multiple shooting
!-------------------------------------------------------------------------
 double precision function dotprod_ms(n_,a,b)
   use variables
   use orbit
   implicit none
   integer,          intent(in) :: n_
   double precision, intent(in) :: a(n*ms), b(n*ms)
   double precision, external :: dotprod
   double precision :: d 
   integer :: p
   
   d = 0d0
   do p = 0, ms-1
      d = d + dotprod(n_,a(p*n+1),b(p*n+1))
   end do
   dotprod_ms = d

 end function dotprod_ms


!-------------------------------------------------------------------------
! Energy norm
!-------------------------------------------------------------------------
 subroutine normE(v, E)
   use variables
   use orbit
   implicit none
   double precision, intent(in)  :: v(n*ms)
   double precision, intent(out) :: E
   double precision :: E_,Ek(0:i_K1),Em(0:i_M1)
   type (coll) :: c(3)
   integer :: i, p
   
   E = 0d0
   do p = 0, ms-1
      call vec2rtz(n,v(p*n+1), c(1),c(2),c(3))
      do i = 1, 3
         call var_coll_norm(c(i), E_,Ek,Em)
         E = E + E_
      end do
   end do
   E = dsqrt(E)
   
 end subroutine normE


!-------------------------------------------------------------------------
!  timestep, shift
!-------------------------------------------------------------------------
 subroutine steporbit(ndts_,x, y)
   use orbit
   use io
   implicit none
   integer,          intent(in)  :: ndts_
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: y(n)
   double precision  :: p(3),st,sz
   integer :: ndt
   
   p = x(:3) * scaleT
   call vec2rtz(n,x, vel_ur,vel_ut,vel_uz)
   if(ndts_/=1) then
      if(Ctnu) then
         if(ctp==1) d_Re = x(4)*ctpsc
         if(ctp==2) d_alpha = x(4)*ctpsc
         call var_precompute()
      end if 
      tim_dt = p(3) / dble(ndts_)
      call vel_matrices()
   end if
   
   do ndt = 1, ndts_
   
      call vel_imposesym()
      call vel_transform()
      call vel_nonlinear()
      if(ndts_/=1 .and. modulo(ndt-1,i_save_rate2)==0) then
         call io_write_energy()
         call io_write_friction()
         call io_write_timestep()
      end if
      call vel_predictor()

      call vel_transform()
      call vel_nonlinear()
      call vel_corrector()
      call vel_imposesym()
      tim_it = 1
      call tim_check_cgce()

      if(ndts_/=1) then
         tim_t    = tim_t    + tim_dt
         tim_step = tim_step + 1
         if(ndt==ndts_) then
            st = -p(1)
            sz = -p(2)
            call var_coll_shift(st,sz,vel_ur, vel_ur)
            call var_coll_shift(st,sz,vel_ut, vel_ut)
            call var_coll_shift(st,sz,vel_uz, vel_uz)
         end if
      end if
      
      if(ndt==ndts_)  &
         call rtz2vec(n,vel_ur,vel_ut,vel_uz, y)

      if(terminate()) then
         call cpu_time(d_stop)
         print*, ' sec/step  = ', (d_stop-d_start)/real(tim_step)
         print*, ' CPU time  = ', int((d_stop-d_start)/6d1), ' mins.'
         call io_closefiles()
         open (99, file='RUNNING')
         close(99, status='delete')
         stop
      end if
      
   end do


 contains
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
!  Termaination conditions
! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
   logical function terminate()
      logical :: file_exist
            
      if(mpi_rnk==0) then
         terminate = .false.
      
         if(tim_step==i_maxtstep) then
            terminate = .true.
            print*, 'maxtstep reached!'
         end if

         if(d_maxt>0d0 .and. tim_t>=d_maxt) then
            terminate = .true.
            print*, 'maxt reached!'
         end if

         call cpu_time(d_stop)
         if(dble(d_stop-d_start)/36d2 >= d_cpuhours) then
            terminate = .true.
            print*, 'cpuhours reached!'
         end if

         if( modulo(tim_step,i_save_rate2)==0) then
            inquire(file='RUNNING', exist=file_exist)
            if(.not. file_exist) then
               terminate = .true.
               print*, 'RUNNING deleted !'
            end if
         end if
      end if

#ifdef _MPI
      call mpi_bcast(terminate,1,mpi_logical, 0,mpi_comm_world,mpi_er)
#endif

   end function terminate

 end subroutine steporbit
 

!-------------------------------------------------------------------------
!  linearised time step for arnoldi
!-------------------------------------------------------------------------
 subroutine multA(x,y0,epsA, y)
   use parameters
   use newton, only : new_x
   use orbit,  only : ndts, n, ms
   implicit none
   double precision, intent(in)  :: x(n), y0(n*ms), epsA
   double precision, intent(out) :: y(n)
   double precision, external :: dotprod
   double precision :: eps
   integer :: p
   
   y = x   
   do p = 0, ms-1
      eps = dsqrt(dotprod(-1,y,y))
      if(eps==0d0)  stop 'multA: eps=0 (1)'
      eps = epsA * dsqrt(dotprod(-1,new_x(p*n+1),new_x(p*n+1))) / eps
      if(eps==0d0)  stop 'multA: eps=0 (2)'
      y = new_x(p*n+1:p*n+n) + eps*y
      y(:8) = new_x(p*n+1:p*n+8)
      call steporbit(ndts(p),y, y)
      y = (y - y0(p*n+1:p*n+n)) / eps
   end do

 end subroutine multA


!-------------------------------------------------------------------------
! get eigenvalues and eigenfunctions about state
!-------------------------------------------------------------------------
#include "./arnoldi.f"
!-------------------------------------------------------------------------
 subroutine getEigen()
   use io
   use orbit 
   use newton
   implicit none
   double precision, external :: dotprod
   double precision, allocatable :: q(:,:)
   double precision :: sv(n), b(n,ncgd*2), wr(ncgd*3), wi(ncgd*3), h(m,m)
   double precision :: y0(n*ms), eps,epsA, st,sz,t, rerr
   double precision :: d1,d2,d3,E0,E,Ek(0:i_K1),Em(0:i_M1)
   integer :: i,j,k,j_, p, ifail

   allocate(q(n,m))
   epsA = epsJ
   rerr = min(rel_err,gtol)
   t    = sum(new_x(3::n))*scaleT
   do p = 0, ms-1
      call steporbit(ndts(p),new_x(p*n+1), y0(p*n+1))
   end do
   						! initial input 
   call vec2rtz(n,new_x, vel_ur,vel_ut,vel_uz)
   st = 0.3d0*(2d0*d_PI)
   sz = 0.3d0*(2d0*d_PI/d_alpha)
   call var_coll_shift(st,sz,vel_ur, vel_ur)
   call var_coll_shift(st,sz,vel_ut, vel_ut)
   call var_coll_shift(st,sz,vel_uz, vel_uz)
   do i = 1, 2
      if(i==1) call var_harmonic(1,0, j,k)
      if(i==2) call var_harmonic(0,1, j,k)
      if(mpi_rnk==k) vel_ur%Re(:,j) = (1d0-mes_D%r(:,2))*mes_D%r(:,2)
      if(mpi_rnk==k) vel_ut%Re(:,j) = (1d0-mes_D%r(:,2))*mes_D%r(:,2)
   end do
   call vel_imposesym()
   call rtz2vec(n,vel_ur,vel_ut,vel_uz, sv)
      						! main arnoldi loop
   k = 0
   do while(.true.)
      sv(:8) = 0d0
      call arnold(n,k,m,ncgd,dotprod,rerr,sv,h,q,b,wr,wi,ifail)
      if(mpi_rnk==0 .and. k>ncgd+2) then
         print*, 'getEigen: k =', k
         do i = 1, ncgd
            print*, 'getEigen: ', real(dlog(dsqrt(wr(i)**2+wi(i)**2))/t)
         end do 
      end if
      if(ifail==-1) then
         if(mpi_rnk==0) print*, ' arnoldi converged!'
         exit
      else if(ifail==0) then
         call multA(sv,y0,epsA, sv)
      else if(ifail==1) then
         print*, 'getEigen: WARNING: arnoldi reached max its'
         exit
      else if(ifail>=2) then
         print*, 'getEigen: arnoldi error'
         deallocate(q)
         return
      end if   
   end do
                                ! magnitude and angle of Floq mult
   wr(ncgd*2+1:) = sqrt(wr(:ncgd)**2+wi(:ncgd)**2)
   wi(ncgd*2+1:) = atan2(wi(:ncgd),wr(:ncgd)) !in (-pi,pi]
                                ! convert Floquet to growth rates
   wr(ncgd+1:ncgd*2) = wr(:ncgd)
   wi(ncgd+1:ncgd*2) = wi(:ncgd)
   call logeig(n,ncgd,b,wr,wi,t,ifail)
   deallocate(q)
      				! save eigvals and Floq multipliers
   if(mpi_rnk==0) then
      open(99,status='unknown',file='arnoldi.dat')
      write(99,*) ' its  = ', k
      write(99,*) ' epsA = ', real(epsA)
      write(99,*) ' T    = ', t
      do i = 1, ncgd
      				! Floq exp, Floq mult, magn & angle
         write(99,'(1I4,6e16.8)') i, wr(i), wi(i),  &
            wr(ncgd+i), wi(ncgd+i), wr(ncgd*2+i), wi(ncgd*2+i)
      end do
      close(99) 
   end if
      				! make norm of eigvecs same as that
                                ! of converged state and save
   tim_t = 0d0
   call vec2rtz(n,new_x,  vel_ur,vel_ut,vel_uz)
   call var_coll_norm(vel_ur, d1,Ek,Em)
   call var_coll_norm(vel_ut, d2,Ek,Em)
   call var_coll_norm(vel_uz, d3,Ek,Em)
   E0 = d1+d2+d3
   io_save1 = 1000
   call io_save_state()
   j = 0
   do i = 1, ncgd
      p = 1
      if(wi(i)/=0d0) p = 2   
      E = 0d0
      do j_ = 1, p
         call vec2rtz(n,b(1,j+j_), vel_ur,vel_ut,vel_uz)
         call var_coll_norm(vel_ur, d1,Ek,Em)
         call var_coll_norm(vel_ut, d2,Ek,Em)
         call var_coll_norm(vel_uz, d3,Ek,Em)
         E = E + d1+d2+d3
      end do
      do j_ = 1, p
         if(E>0d0)  b(:,j+j_) = b(:,j+j_) * dsqrt(E0/E)
         call vec2rtz(n,b(1,j+j_), vel_ur,vel_ut,vel_uz)
         io_save1 = 1000 + j+j_
         call io_save_state()
      end do
      j = j + p
   end do

 end subroutine getEigen 
