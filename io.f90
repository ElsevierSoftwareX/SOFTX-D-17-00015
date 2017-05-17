!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module io
!**************************************************************************
   use velocity
   use netcdf
   implicit none
   save

   character(200)        :: io_statefile
   integer               :: io_save1, io_save2
   integer,     private  :: io_KE, io_ID, io_dt, io_pt, io_fr
   type (coll), private  :: c1,c2,c3

 contains
 
!--------------------------------------------------------------------------
!  initialiser fn
!--------------------------------------------------------------------------
   subroutine io_precompute()
      io_statefile = 'state.cdf.in'
      io_save1 = 0
      io_save2 = 0
      io_dt    = 10
      io_KE    = 20
      io_ID    = 21
      io_pt    = 0
      io_fr    = 50
   end subroutine io_precompute 
 

!--------------------------------------------------------------------------
!  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   subroutine io_openfiles()
      character(10), save :: s = 'unknown', a = 'sequential'
      if(mpi_rnk/=0) return
      if(io_dt/=0)  open(io_dt,status=s,access=a, file='tim_step.dat')
      if(io_KE/=0)  open(io_KE,status=s,access=a, file='vel_energy.dat')
      if(io_ID/=0)  open(io_ID,status=s,access=a, file='vel_totEID.dat')
      if(io_pt/=0)  open(io_pt,status=s,access=a, file='vel_point.dat')
      if(io_fr/=0)  open(io_fr,status=s,access=a, file='vel_friction.dat')
!      s = 'old'
      a = 'append'
   end subroutine io_openfiles


!--------------------------------------------------------------------------
!  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      if(mpi_rnk/=0) return
      if(io_dt/=0) close(io_dt)
      if(io_KE/=0) close(io_KE)
      if(io_ID/=0) close(io_ID)
      if(io_pt/=0) close(io_pt)
      if(io_fr/=0) close(io_fr)
   end subroutine io_closefiles


!--------------------------------------------------------------------------
!  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files()

      if(modulo(tim_step,i_save_rate1)==0) then
         call io_save_state()
         call io_save_spectrum()
         call io_save_meanprof()
         io_save1 = io_save1+1
      endif

      if(modulo(tim_step,i_save_rate2)==0) then
         if(io_KE/=0) call io_write_energy()
         if(io_ID/=0) call io_write_totEID()
         if(io_pt/=0) call io_write_pointvel()
         if(io_fr/=0) call io_write_friction()
         if(io_dt/=0 .and. d_timestep>0d0) call io_write_timestep()
         io_save2 = io_save2+1
      end if

      if(io_dt/=0 .and. tim_new_dt) call io_write_timestep()

      if(modulo(tim_step,i_save_rate2*50)==0) then
         call io_closefiles()
         call io_openfiles()
      end if

   end subroutine io_write2files


!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state()
      integer :: e, f, i, rd
      integer :: n, N_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: A(lda,lda,i_N)
      double precision, allocatable :: r(:)
      logical :: interp

      e=nf90_open(io_statefile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
         if(mpi_rnk==0) print*, 'state file not found!: '//io_statefile
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'io_load_state: file not found!'
      end if

      e=nf90_get_att(f,nf90_global,'t', d)
      if(d_time<0d0) tim_t = d
      if(mpi_rnk==0 .and. dabs(tim_t-d)>1d-8)  &
         print*,' t    :',d,' --> ', tim_t

      e=nf90_get_att(f,nf90_global,'Re', d)
      if(mpi_rnk==0 .and. dabs(d_Re-d)>1d-8)  &
         print*,' Re   :',d,' --> ', d_Re
      e=nf90_get_att(f,nf90_global,'alpha', d)
      if(mpi_rnk==0 .and. dabs(d_alpha-d)>1d-8)  &
         print*,' alpha:',d,' --> ', d_alpha

      e=nf90_inq_dimid(f,'r', rd)
      e=nf90_inquire_dimension(f,rd, len=N_)

      allocate( r(1-i_KL:N_))
      e=nf90_inq_varid(f,'r', i)
      e=nf90_get_var(f,i, r(1:N_))
      r(1-i_KL:0) = -r(i_KL:1:-1)
      interp = (N_/=i_N)
      if(.not.interp) interp = (maxval(dabs(mes_D%r(:,1)-r(1:N_)))>1d-8)
      if(interp) then
         if(mpi_rnk==0) print*,' N    :', N_, ' --> ',i_N 
         call io_interp_wts(i_KL+N_,r,i_N,mes_D%r(1,1), A)
      end if

      e=nf90_inq_varid(f,'dt', i)
      e=nf90_get_var(f,i, d)
      if(d_timestep<0d0) tim_dt = d
      if(d_timestep>0d0) tim_dt = d_timestep
      e=nf90_inq_varid(f,'dtcor', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_corr_dt)
      e=nf90_inq_varid(f,'dtcfl', i)
      if(e==nf90_noerr)  e=nf90_get_var(f,i, tim_cfl_dt)

      call io_load_coll(f,'Ur',interp,N_,r,A,1, vel_ur)
      call io_load_coll(f,'Ut',interp,N_,r,A,1, vel_ut)
      call io_load_coll(f,'Uz',interp,N_,r,A,0, vel_uz)      

      deallocate(r)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_state


!--------------------------------------------------------------------------
!  Load coll variable
!--------------------------------------------------------------------------
   subroutine io_load_coll(f,nm,interp,N__,r,W,S, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      logical,          intent(in)  :: interp
      integer,          intent(in)  :: N__
      double precision, intent(in)  :: r(i_KL+N__)
      double precision, intent(in)  :: W(2*i_KL+1,2*i_KL+1,*)
      integer,          intent(in)  :: S
      type (coll),      intent(out) :: a
      double precision :: fn(1-i_KL:N__)
      integer :: K1,M1, nh,nh_,nh__,n,k,m
      integer :: K__, M__, Mp__      
      integer :: e,i
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'M',  M__)
      e=nf90_get_att(f,i,'Mp', Mp__)
      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_M)  print*, nm, ' M :', M__, ' --> ',i_M 
         if(Mp__/=i_Mp) print*, nm, ' Mp:', Mp__,' --> ',i_Mp 
      end if
      if(Mp__/=i_Mp .and. M__>2) stop 'io_load_coll: Mp /= i_Mp'

      a%Re = 0d0
      a%Im = 0d0

      K1 = max(K__,i_K)-1
      M1 = max(M__,i_M)-1

      nh__ = -1
      nh = -1
      do m = 0, M1
         do k = -K1, K1
            if(k<0 .and. m==0) cycle
            if(abs(k)>=i_K .or. m>=i_M)  nh__ = nh__ + 1
            if(abs(k)>=i_K .or. m>=i_M)  cycle
            nh   = nh + 1
            if(abs(k)>=K__ .or. m>=M__)  cycle
            nh__ = nh__ + 1
            if(nh<var_H%pH0 .or. nh>var_H%pH0+var_H%pH1)  cycle
            nh_  = nh-var_H%pH0
            if(interp) then
               e=nf90_get_var(f,i, fn(1:N__),start=(/1,nh__+1,1/))
               fn(:0) = fn(i_KL:1:-1)
               if(modulo(m*i_Mp+S,2)==1)  fn(:0) = -fn(:0)
               call io_interp(i_KL+N__,r,fn,W,i_N,mes_D%r(1,1), a%Re(1,nh_))
               e=nf90_get_var(f,i, fn(1:N__),start=(/1,nh__+1,2/))
               fn(:0) = fn(i_KL:1:-1)
               if(modulo(m*i_Mp+S,2)==1)  fn(:0) = -fn(:0)
               call io_interp(i_KL+N__,r,fn,W,i_N,mes_D%r(1,1), a%Im(1,nh_))
            else
               e=nf90_get_var(f,i, a%Re(1:N__,nh_),start=(/1,nh__+1,1/))
               e=nf90_get_var(f,i, a%Im(1:N__,nh_),start=(/1,nh__+1,2/))
            end if
         end do
      end do

   end subroutine io_load_coll


!--------------------------------------------------------------------------
!  interpolate; return weights
!--------------------------------------------------------------------------
   subroutine io_interp_wts(ni,xi,no,xo, A)
      integer, parameter :: lda = 2*i_KL+1
      integer,          intent(in)  :: ni, no
      double precision, intent(in)  :: xi(ni), xo(no)
      double precision, intent(out) :: A(lda,lda,no) 
      integer :: n,nn,i,j,l,r
      
      do n = 1, no
         j = 1
         do while(xi(j)<xo(n)-1d-8 .and. j<ni)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         do i = 1, nn
            A(i,1,n) = 1d0
         end do
         do j = 2, nn
            do i = 1, nn
               A(i,j,n) = A(i,j-1,n) * (xi(l+i-1)-xo(n)) / dble(j-1) 
            end do
         end do
         call mes_mat_invert(nn,A(1,1,n),lda)
      end do

   end subroutine io_interp_wts


!--------------------------------------------------------------------------
!  interpolate, given weights from io_interp_wts()
!--------------------------------------------------------------------------
   subroutine io_interp(ni,xi,fi,A,no,xo, fo)
      integer, parameter :: lda = 2*i_KL+1
      integer,          intent(in)  :: ni, no
      double precision, intent(in)  :: xi(ni), fi(ni), xo(no)
      double precision, intent(in)  :: A(lda,lda,no) 
      double precision, intent(out) :: fo(no)
      integer :: n,nn,i,j,l,r
      
      do n = 1, no
         j = 1
         do while(xi(j)<xo(n)-1d-8 .and. j<ni)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         fo(n) = dot_product(A(1,1:nn,n),fi(l:r))
      end do

   end subroutine io_interp


!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state()
      character(4) :: cnum
      integer :: e, f
      integer :: rd, Hd, ReImd, dims(3)
      integer :: r,dt,dtcor,dtcfl, Ur,Ut,Uz

      write(cnum,'(I4.4)') io_save1

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t
         e=nf90_create('state'//cnum//'.cdf.dat', nf90_clobber, f)

         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)

         e=nf90_def_dim(f, 'r', i_N, rd)
         e=nf90_def_dim(f, 'H', i_H1+1, Hd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         e=nf90_def_var(f, 'r',     nf90_double, (/rd/), r)
         e=nf90_def_var(f, 'dt',    nf90_double, dt)
         e=nf90_def_var(f, 'dtcor', nf90_double, dtcor)
         e=nf90_def_var(f, 'dtcfl', nf90_double, dtcfl)

         dims = (/rd,Hd,ReImd/)
         call io_define_coll(f, 'Ur', dims, Ur)
         call io_define_coll(f, 'Ut', dims, Ut)         
         call io_define_coll(f, 'Uz', dims, Uz)         

         e=nf90_enddef(f)

         e=nf90_put_var(f, r, mes_D%r(1:i_N,1))
         e=nf90_put_var(f, dt, tim_dt)
         e=nf90_put_var(f, dtcor, tim_corr_dt)
         e=nf90_put_var(f, dtcfl, tim_cfl_dt)
      end if

      call io_save_coll(f,Ur, vel_ur)
      call io_save_coll(f,Ut, vel_ut)
      call io_save_coll(f,Uz, vel_uz)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_state
 

!--------------------------------------------------------------------------
!  Save coll variable
!--------------------------------------------------------------------------
   subroutine io_define_coll(f,nm,dims, id)
      integer,      intent(in) :: f, dims(3)
      character(*), intent(in) :: nm
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'K', i_K)      
      e=nf90_put_att(f, id,  'M', i_M)
      e=nf90_put_att(f, id, 'Mp', i_Mp)
   end subroutine io_define_coll

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine io_save_coll(f,id,a)
      integer,     intent(in) :: f, id
      type (coll), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_N,0:i_H1), start=(/1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_N,0:i_H1), start=(/1,1,2/))

#else
      integer :: r, pH0,pH1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_N,0:var_H%pH1), start=(/1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_N,0:var_H%pH1), start=(/1,1,2/))         
         do r = 1, mpi_sze-1
            pH0 = var_H%pH0_(r)
            pH1 = var_H%pH1_(r)
            mpi_tg = r
            call mpi_recv( c1%Re, i_N*(pH1+1), mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( c1%Im, i_N*(pH1+1), mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,c1%Re(1:i_N,0:pH1), start=(/1,pH0+1,1/))
            e=nf90_put_var(f,id,c1%Im(1:i_N,0:pH1), start=(/1,pH0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re, i_N*(var_H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im, i_N*(var_H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_coll


!--------------------------------------------------------------------------
!  save spectrum
!--------------------------------------------------------------------------
   subroutine io_save_spectrum()
      double precision :: n_(1:i_N), k_(0:i_K1), m_(0:i_M1)
      double precision :: n__(1:i_N), k__(0:i_K1), m__(0:i_M1)
      double precision,save :: TM(i_N,i_N), x(i_N)
      double precision :: d(i_N), dRe(i_N), dIm(i_N)
      logical, save :: set=.false.
      character(4) :: cnum
      integer :: i, n,kp
      _loop_km_vars
   10 format(i4,1e20.12)
      
      if(.not.set) then
         set =.true.
         do n = 0, i_N-1
            x(n+1) = 0.5d0 * ( 1d0 + dcos(d_PI*(i_N-n)/dble(i_N)) )
         end do
         do n = 1, i_N
            call cheby(n-1, 0, x, i_N, TM(1,n))
         end do
         call mes_mat_invert(i_N,TM,i_N)
         TM = transpose(TM)
      end if

      n_ = 0d0
      k_ = 0d0
      m_ = 0d0
      _loop_km_begin
         dRe = matmul(vel_ur%Re(:,nh), TM)
         dIm = matmul(vel_ur%Im(:,nh), TM)
         d = dRe*dRe+dIm*dIm
         dRe = matmul(vel_ut%Re(:,nh), TM)
         dIm = matmul(vel_ut%Im(:,nh), TM)
         d = max(d, dRe*dRe+dIm*dIm)
         dRe = matmul(vel_uz%Re(:,nh), TM)
         dIm = matmul(vel_uz%Im(:,nh), TM)
         d = max(d, dRe*dRe+dIm*dIm)
         d = dsqrt(d)
         kp = abs(k)
         do n = 1, i_N
             n_(n)  = max(d(n), n_(n))
             k_(kp) = max(d(n), k_(kp))
             m_(m)  = max(d(n), m_(m))            
         end do
      _loop_km_end
#ifdef _MPI
      call mpi_allreduce(n_, n__, i_N, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      n_ = n__
      call mpi_allreduce(k_, k__, i_K, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      k_ = k__
      call mpi_allreduce(m_, m__, i_M, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      m_ = m__
#endif

      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_spec'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# n'
      do i = 1, i_N
         write(11,10) i, n_(i)      
      end do
      write(11,*)
      write(11,*) '# k'
      do i = 0, i_K1
         write(11,10) i, k_(i)      
      end do
      write(11,*)
      write(11,*) '# m'
      do i = 0, i_M1
         write(11,10) i*i_Mp, m_(i)      
      end do
      close(11)

   end subroutine io_save_spectrum


!--------------------------------------------------------------------------
!  save axial mean flow profile
!--------------------------------------------------------------------------
   subroutine io_save_meanprof()
      character(4) :: cnum
      integer :: n

      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_prof'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# r  uz(r)'
      do n = 1, i_N
         write(11,'(2e20.12)')  mes_D%r(n,1),  &
            vel_uz%Re(n,0) + 1d0-mes_D%r(n,2)      
      end do
      close(11)

   end subroutine io_save_meanprof

 
!--------------------------------------------------------------------------
!  write to energy file
!--------------------------------------------------------------------------
   subroutine io_write_energy()
      double precision :: E,Ek0,Em0,E_, Ek(0:i_K1), Em(0:i_M1)

      call var_coll_norm(vel_ur, E,Ek,Em)
      Ek0 = Ek(0)
      Em0 = Em(0)
      call var_coll_norm(vel_ut, E_,Ek,Em)
      E   = E + E_
      Ek0 = Ek0 + Ek(0)
      Em0 = Em0 + Em(0)
      call var_coll_norm(vel_uz, E_,Ek,Em)
      E   = E + E_
      Ek0 = Ek0 + Ek(0)
      Em0 = Em0 + Em(0)
      
      if(mpi_rnk/=0) return
      write(io_KE,'(4e20.12)')  tim_t, E, Ek0, Em0
      
      if(E-Ek0>d_minE3d .or. tim_t<20d0) return
      print*, 'io_write_energy: Relaminarised!'
      open(99,file='RUNNING')
      close(99, status='delete')

   end subroutine io_write_energy


!--------------------------------------------------------------------------
!  (Total) Energy/E_lam, Input/I_lam, Dissip/D_lam
!--------------------------------------------------------------------------
   subroutine io_write_totEID()
      double precision :: E,E_,E__, Ek(0:i_K1), Em(0:i_M1)
      double precision :: Inp, Diss, Enrg, Lz

      call var_coll_copy(vel_uz, c3)
      if(mpi_rnk==0) c3%Re(:,0) = c3%Re(:,0) + vel_U
      call var_coll_norm(vel_ur, E,Ek,Em)
      call var_coll_norm(vel_ut, E_,Ek,Em)
      call var_coll_norm(c3,     E__,Ek,Em)
      Enrg = E + E_ + E__
      
      call var_coll_curl(vel_ur,vel_ut,c3, c1,c2,c3)
      call var_coll_norm(c1, E,Ek,Em)
      call var_coll_norm(c2, E_,Ek,Em)
      call var_coll_norm(c3, E__,Ek,Em)
      Diss = E + E_ + E__
      
      Lz = 2d0*d_PI/d_alpha
      Enrg = Enrg / (d_PI**2/(3d0*d_alpha))
      Diss = Diss * 2d0/d_Re
      Diss = Diss / abs(vel_Up(i_N)*d_PI*Lz/d_Re)
                                                   !   I/I_lam = 1+beta
      if(mpi_rnk/=0) return
      Inp = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      Inp = Inp + vel_Up(i_N)
      Inp = abs(Inp/vel_Up(i_N))
      
      write(io_ID,'(4e20.12)')  tim_t, Enrg, Inp, Diss

   end subroutine io_write_totEID


!--------------------------------------------------------------------------
!  write velocity at points.  r,t,z components at 3 points.
!--------------------------------------------------------------------------
   subroutine io_write_pointvel()
      integer :: n,rt(3),r(3)
      double precision :: d(9), d_
      if(_Ns/=1) stop 'io_write_pointvel: put _Ns=1'

      d = 1d0      
      do n = 1, i_N
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.6667d0)
         if(d_<d(1)) r(1) = n 
         if(d_<d(1)) d(1) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.7550d0)
         if(d_<d(2)) r(2) = n 
         if(d_<d(2)) d(2) = d_
         d_ = dabs(mes_D%r(n,1)-0.6d0) !0.9217d0)
         if(d_<d(3)) r(3) = n 
         if(d_<d(3)) d(3) = d_
      end do
      do n = 0, mpi_sze-1
         if(mes_D%pNi_(n)<=r(1)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(1)) rt(1) = n
         if(mes_D%pNi_(n)<=r(2)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(2)) rt(2) = n
         if(mes_D%pNi_(n)<=r(3)  &
            .and. mes_D%pNi_(n)+mes_D%pN_(n)-1>=r(3)) rt(3) = n
      end do

      if(rt(1)==mpi_rnk) then
         d(1) = vel_r%Re(0,0,r(1)-mes_D%pNi+1)
         d(2) = vel_t%Re(0,0,r(1)-mes_D%pNi+1)
         d(3) = vel_z%Re(0,0,r(1)-mes_D%pNi+1)
      end if      
      if(rt(2)==mpi_rnk) then
         d(4) = vel_r%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(5) = vel_t%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
         d(6) = vel_z%Re(i_pZ/3,i_Th/3,r(2)-mes_D%pNi+1)
      end if
      if(rt(3)==mpi_rnk) then
         d(7) = vel_r%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(8) = vel_t%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
         d(9) = vel_z%Re(2*i_pZ/3,2*i_Th/3,r(3)-mes_D%pNi+1)
      end if

#ifdef _MPI
      call mpi_bcast(d(1), 3, mpi_double_precision,  &
         rt(1), mpi_comm_world, mpi_er)
      call mpi_bcast(d(4), 3, mpi_double_precision,  &
         rt(2), mpi_comm_world, mpi_er)
      call mpi_bcast(d(7), 3, mpi_double_precision,  &
         rt(3), mpi_comm_world, mpi_er)
#endif
      if(mpi_rnk/=0) return
      if(tim_step==0)  write(io_pt,*) '# r=,', real(mes_D%r(r(1),1)),  &
         real(mes_D%r(r(2),1)), real(mes_D%r(r(3),1))
      write(io_pt,'(10e16.8)') tim_t, (d(n),n=1,9)

   end subroutine io_write_pointvel


!--------------------------------------------------------------------------
!  write:  1, time;  2,  bulk vel / excess pressure fraction if fixed flux;
!          3, mean uz at r=0;  4, friction velocity. 
!--------------------------------------------------------------------------
   subroutine io_write_friction()
      double precision :: Ub, Uc, Ufr

      Ub = 0.5d0 + 2d0*dot_product(vel_uz%Re(:,0),mes_D%intrdr)
      if(b_const_flux)  Ub = vel_Pr0
      Uc = 1d0 + dot_product(vel_uz%Re(1:1+i_KL,0),mes_D%dr0(:,0))
      Ufr = dot_product(vel_uz%Re(i_N-i_KL:i_N,0),mes_D%dr1(:,1))
      Ufr = dsqrt(dabs((Ufr-2d0)/d_Re))
      
      if(mpi_rnk/=0) return
      write(io_fr,'(4e20.12)') tim_t, Ub, Uc, Ufr
   
   end subroutine io_write_friction


!--------------------------------------------------------------------------
!  write to timestep file
!--------------------------------------------------------------------------
   subroutine io_write_timestep()
   10 format(3e18.10,I2)
   11 format(1e18.10,I8,3e17.9,I2)
      if(tim_new_dt) then
         if(mpi_rnk/=0) return
         write(io_dt,11) tim_t, tim_step,  &
            tim_dt, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      else
         call vel_maxtstep()
         if(mpi_rnk/=0) return
         write(io_dt,10) tim_t, tim_corr_dt, tim_cfl_dt, tim_cfl_dir
      end if
   end subroutine io_write_timestep


!**************************************************************************
 end module io
!**************************************************************************

