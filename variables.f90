!*************************************************************************
!  Variables in spectral,collocated and physical space.
!  Simple functions on variables.
!
!     r_n    n in [1,N]
!     th_m   m in [0,Th)
!     z_k    k in [0,Z)
!
!  harmonics ordered s.t. usual loop:
!     nh = -1
!     do m = 0, M-1             nh  harmonic index
!        m_ = m*Mp              m   m index
!        do k = -(K-1), K-1     m_  actual m mode
!           if(m==0 .and. k<0)  cycle
!           nh = nh + 1         k   k mode
!           ...
!
!*************************************************************************
#include "../parallel.h"
 module variables
!*************************************************************************
   use mpif
   use parameters
   use meshs
   implicit none
   save

   type harm
      integer              :: pH0,pH1, pH0_(0:_Np-1),pH1_(0:_Np-1)
   end type harm

   type spec
      double precision     :: Re(0:_Hs1, i_pN)
      double precision     :: Im(0:_Hs1, i_pN)
   end type spec

   type coll
      double precision     :: Re(i_N, 0:i_pH1)
      double precision     :: Im(i_N, 0:i_pH1)
   end type coll

   type phys
      double precision     :: Re(0:i_pZ-1, 0:i_Th-1, i_pN)
   end type phys
   

   type (harm)               :: var_H
   type (coll),      private :: c1,c2,c3
   double precision, private :: ad_k2a2(-i_K1:i_K1)
   double precision, private :: ad_k1a1(-i_K1:i_K1)
   double precision, private :: ad_m1r1(i_N,0:i_M1)
   double precision, private :: ad_m2r2(i_N,0:i_M1)

 contains


!------------------------------------------------------------------------
!  Initialise stuff for variables
!------------------------------------------------------------------------
   subroutine var_precompute()
      integer :: n,r, ms,H1
         			! distribute modes accross processors
      var_H%pH1_ = -1
      do ms = 0, _Ns-1
         H1 = _Hs1
         if(ms==0) H1 = H1-i_K1
         do n = 0, H1
            r = _Nr - modulo(n,_Nr) - 1 + ms*_Nr
            var_H%pH1_(r) = var_H%pH1_(r) + 1
         end do
      end do
      var_H%pH0_(0) = 0
      do r = 1, _Np-1
         var_H%pH0_(r) = var_H%pH0_(r-1) + var_H%pH1_(r-1) + 1
      end do
      var_H%pH0 = var_H%pH0_(mpi_rnk)
      var_H%pH1 = var_H%pH1_(mpi_rnk)
      
      do n = 0, i_M1
         ad_m1r1(:,n) = mes_D%r(:,-1) * i_Mp*n
         ad_m2r2(:,n) = mes_D%r(:,-2) * i_Mp*n * i_Mp*n
      end do
      do n = -i_K1, i_K1
         ad_k2a2(n) = d_alpha*n * d_alpha*n
         ad_k1a1(n) = d_alpha*n
      end do

   end subroutine var_precompute


!-------------------------------------------------------------------------
!  get nh and mpi_rank for given k,m
!-------------------------------------------------------------------------
   subroutine var_harmonic(k,m, nh,p)
      integer, intent(in)  :: k,m
      integer, intent(out) :: nh,p
      nh = m*(2*i_K1+1) + k
      do p = (m/_Ms)*_Nr, mpi_sze-1
         if(nh>=var_H%pH0_(p) .and. nh<=var_H%pH0_(p)+var_H%pH1_(p)) then
            nh = nh-var_H%pH0_(p)
            return
         end if
      end do
   end subroutine var_harmonic


!-------------------------------------------------------------------------
!  null function
!-------------------------------------------------------------------------
   subroutine var_null(flag)
      integer, intent(in), optional :: flag
   end subroutine var_null


!-------------------------------------------------------------------------
!  initialise collocated variable
!-------------------------------------------------------------------------
   subroutine var_coll_init(a)
      type (coll), intent(out) :: a
      a%Re = 0d0
      a%Im = 0d0
   end subroutine var_coll_init


!-------------------------------------------------------------------------
!  Copy a collocated variable
!-------------------------------------------------------------------------
   subroutine var_coll_copy(in, out)
      type (coll), intent(in)  :: in
      type (coll), intent(out) :: out
      out%Re(:,0:var_H%pH1) = in%Re(:,0:var_H%pH1)
      out%Im(:,0:var_H%pH1) = in%Im(:,0:var_H%pH1)
   end subroutine var_coll_copy


!------------------------------------------------------------------------
!     out := out + in
!------------------------------------------------------------------------
   subroutine var_coll_add(ac, a)
      type (coll), intent(in)    :: ac
      type (coll), intent(inout) :: a
      a%Re(:,0:var_H%pH1) = a%Re(:,0:var_H%pH1) + ac%Re(:,0:var_H%pH1)
      a%Im(:,0:var_H%pH1) = a%Im(:,0:var_H%pH1) + ac%Im(:,0:var_H%pH1)
   end subroutine var_coll_add


!------------------------------------------------------------------------
!     out := out - in
!------------------------------------------------------------------------
   subroutine var_coll_sub(ac, a)
      type (coll), intent(in)    :: ac
      type (coll), intent(inout) :: a
      a%Re(:,0:var_H%pH1) = a%Re(:,0:var_H%pH1) - ac%Re(:,0:var_H%pH1)
      a%Im(:,0:var_H%pH1) = a%Im(:,0:var_H%pH1) - ac%Im(:,0:var_H%pH1)
   end subroutine var_coll_sub


!------------------------------------------------------------------------
!     out := d * in
!------------------------------------------------------------------------
   subroutine var_coll_dmult(d, a)
      double precision, intent(in) :: d
      type (coll),   intent(inout) :: a
      a%Re(:,0:var_H%pH1) = a%Re(:,0:var_H%pH1) * d
      a%Im(:,0:var_H%pH1) = a%Im(:,0:var_H%pH1) * d
   end subroutine var_coll_dmult


!------------------------------------------------------------------------
!     out := - in
!------------------------------------------------------------------------
   subroutine var_coll_neg(ac, a)
      type (coll), intent(in)    :: ac
      type (coll), intent(inout) :: a
      a%Re(:,0:var_H%pH1) = - ac%Re(:,0:var_H%pH1)
      a%Im(:,0:var_H%pH1) = - ac%Im(:,0:var_H%pH1)
   end subroutine var_coll_neg


!-------------------------------------------------------------------------
!  convert collocated -> spectral
!-------------------------------------------------------------------------
#ifndef _MPI
   subroutine var_coll2spec(c,s, c2,s2, c3,s3)
      type (coll), intent(in)  :: c
      type (spec), intent(out) :: s
      type (coll), intent(in),  optional :: c2,c3
      type (spec), intent(out), optional :: s2,s3
      s%Re = transpose(c%Re)
      s%Im = transpose(c%Im)
      if(.not. present(c2)) return
      s2%Re = transpose(c2%Re)
      s2%Im = transpose(c2%Im)
      if(.not. present(c3)) return
      s3%Re = transpose(c3%Re)
      s3%Im = transpose(c3%Im)
   end subroutine var_coll2spec

#else
   subroutine var_coll2spec(c,s, c2,s2, c3,s3)
      type (coll), intent(in)  :: c
      type (spec), intent(out) :: s
      type (coll), intent(in),  optional :: c2,c3
      type (spec), intent(out), optional :: s2,s3
      double precision :: bsend(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
      double precision :: brecv(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
      integer :: stp, dst,src, n,nh,l, nc, rko,nho

      nc = 1
      if(present(c2)) nc = 2
      if(present(c3)) nc = 3
  
      rko = (mpi_rnk/_Nr)*_Nr
      nho = var_H%pH0_(rko)

      do stp = 0, _Nr-1
         src  = modulo(mpi_sze-stp+mpi_rnk, _Nr) + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*mes_D%pN*(var_H%pH1_(src)+1)*nc,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nr-1
         dst  = modulo(stp+mpi_rnk, _Nr) + rko
         l = 1
         do n = mes_D%pNi_(dst), mes_D%pNi_(dst)+mes_D%pN_(dst)-1
            do nh = 0, var_H%pH1
               bsend(l,  stp) = c%Re(n,nh)
               bsend(l+1,stp) = c%Im(n,nh)
               l = l + 2
               if(nc<2) cycle
               bsend(l,  stp) = c2%Re(n,nh)
               bsend(l+1,stp) = c2%Im(n,nh)
               l = l + 2
               if(nc<3) cycle
               bsend(l,  stp) = c3%Re(n,nh)
               bsend(l+1,stp) = c3%Im(n,nh)
               l = l + 2
            end do
         end do
         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*mes_D%pN_(dst)*(var_H%pH1+1)*nc,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do
      
      do stp = 0, _Nr-1
         src  = modulo(mpi_sze-stp+mpi_rnk, _Nr) + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do n = 1, mes_D%pN
            do nh =  &
                  var_H%pH0_(src)-nho, var_H%pH0_(src)+var_H%pH1_(src)-nho
               s%Re(nh,n) = brecv(l,  stp)
               s%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
               if(nc<2) cycle
               s2%Re(nh,n) = brecv(l,  stp)
               s2%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
               if(nc<3) cycle
               s3%Re(nh,n) = brecv(l,  stp)
               s3%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
            end do
         end do
      end do

      do stp = 0, _Nr-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do

   end subroutine var_coll2spec
#endif
   
!-------------------------------------------------------------------------
!  convert spectral -> collocated
!-------------------------------------------------------------------------
#ifndef _MPI
   subroutine var_spec2coll(s,c, s2,c2, s3,c3)
      type (spec), intent(in)  :: s
      type (coll), intent(out) :: c
      type (spec), intent(in),  optional :: s2,s3
      type (coll), intent(out), optional :: c2,c3
      c%Re = transpose(s%Re)
      c%Im = transpose(s%Im)
      if(.not. present(s2)) return
      c2%Re = transpose(s2%Re)
      c2%Im = transpose(s2%Im)
      if(.not. present(s3)) return
      c3%Re = transpose(s3%Re)
      c3%Im = transpose(s3%Im)
   end subroutine var_spec2coll

#else
   subroutine var_spec2coll(s,c, s2,c2, s3,c3)
      type (spec), intent(in)  :: s
      type (coll), intent(out) :: c
      type (spec), intent(in),  optional :: s2,s3
      type (coll), intent(out), optional :: c2,c3
      double precision :: bsend(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
      double precision :: brecv(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
      integer :: stp, dst,src, n,nh,l, ns, rko,nho

      ns = 1
      if(present(s2)) ns = 2
      if(present(s3)) ns = 3

      rko = (mpi_rnk/_Nr)*_Nr
      nho = var_H%pH0_(rko)

      do stp = 0, _Nr-1
         src  = modulo(mpi_sze-stp+mpi_rnk, _Nr) + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*(var_H%pH1+1)*mes_D%pN_(src)*ns,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nr-1
         dst  = modulo(stp+mpi_rnk, _Nr) + rko
         l = 1
         do nh = var_H%pH0_(dst)-nho, var_H%pH0_(dst)+var_H%pH1_(dst)-nho
            do n = 1, mes_D%pN
               bsend(l,  stp) = s%Re(nh,n)
               bsend(l+1,stp) = s%Im(nh,n)
               l = l + 2
               if(ns<2) cycle
               bsend(l,  stp) = s2%Re(nh,n)
               bsend(l+1,stp) = s2%Im(nh,n)
               l = l + 2
               if(ns<3) cycle
               bsend(l,  stp) = s3%Re(nh,n)
               bsend(l+1,stp) = s3%Im(nh,n)
               l = l + 2
            end do
         end do
         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*(var_H%pH1_(dst)+1)*mes_D%pN*ns,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do
      
      do stp = 0, _Nr-1
         src  = modulo(mpi_sze-stp+mpi_rnk, _Nr) + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do nh = 0, var_H%pH1
            do n = mes_D%pNi_(src), mes_D%pNi_(src)+mes_D%pN_(src)-1
               c%Re(n,nh) = brecv(l,  stp)
               c%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
               if(ns<2) cycle
               c2%Re(n,nh) = brecv(l,  stp)
               c2%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
               if(ns<3) cycle
               c3%Re(n,nh) = brecv(l,  stp)
               c3%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
            end do
         end do
      end do

      do stp = 0, _Nr-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do

   end subroutine var_spec2coll
#endif


!------------------------------------------------------------------------
!  Multiply a collocated variable by a mesh-type
!     out = A.in
!------------------------------------------------------------------------
   subroutine var_coll_meshmult(S,A,in, out)
      integer,     intent(in)  :: S
      type (mesh), intent(in)  :: A
      type (coll), intent(in)  :: in
      type (coll), intent(out) :: out
      double precision :: re(i_N), im(i_N)
      integer :: n,j,l,r
      _loop_km_vars

      _loop_km_begin
         re = 0d0
         im = 0d0
         do j = 1, i_N
            l = max(1,j-i_KL)
            r = min(j+i_KL,i_N)
            do n = l, r
               re(n) = re(n) + A%M(i_KL+1+n-j,j) * in%Re(j,nh)
               im(n) = im(n) + A%M(i_KL+1+n-j,j) * in%Im(j,nh)
            end do
         end do
         out%Re(:,nh) = re
         out%Im(:,nh) = im
      _loop_km_end

   end subroutine var_coll_meshmult


!------------------------------------------------------------------------
!  take the curl of a vector
!------------------------------------------------------------------------
   subroutine var_coll_curl(r,t,z, or,ot,oz)
      type (coll), intent(in)  :: r,t,z
      type (coll), intent(out) :: or,ot,oz
      _loop_km_vars

      call var_coll_copy(r,c1)
      call var_coll_copy(t,c2)
      call var_coll_copy(z,c3)
 
      _loop_km_begin
         or%Re(:,nh) = -c3%Im(:,nh)*ad_m1r1(:,m) + c2%Im(:,nh)*ad_k1a1(k)
         or%Im(:,nh) =  c3%Re(:,nh)*ad_m1r1(:,m) - c2%Re(:,nh)*ad_k1a1(k)
      _loop_km_end

      call var_coll_meshmult(0,mes_D%dr(1),c3, ot)
      _loop_km_begin
         ot%Re(:,nh) = -c1%Im(:,nh)*ad_k1a1(k) - ot%Re(:,nh)
         ot%Im(:,nh) =  c1%Re(:,nh)*ad_k1a1(k) - ot%Im(:,nh)
      _loop_km_end

      call var_coll_meshmult(1,mes_D%dr(1),c2, oz)
      _loop_km_begin
         oz%Re(:,nh) =  c2%Re(:,nh)*mes_D%r(:,-1) + oz%Re(:,nh)  &
                       +c1%Im(:,nh)*ad_m1r1(:,m)
         oz%Im(:,nh) =  c2%Im(:,nh)*mes_D%r(:,-1) + oz%Im(:,nh)  &
                       -c1%Re(:,nh)*ad_m1r1(:,m)
      _loop_km_end

   end subroutine var_coll_curl


!------------------------------------------------------------------------
!  take the gradient of a scalar
!------------------------------------------------------------------------
   subroutine var_coll_grad(p, r,t,z)
      type (coll), intent(in)  :: p
      type (coll), intent(out) :: r,t,z
      _loop_km_vars

      call var_coll_copy(p,c1)

      call var_coll_meshmult(1,mes_D%dr(1),c1, r)
      _loop_km_begin
         t%Re(:,nh) = -c1%Im(:,nh)*ad_m1r1(:,m)
         t%Im(:,nh) =  c1%Re(:,nh)*ad_m1r1(:,m)
      _loop_km_end
      _loop_km_begin
         z%Re(:,nh) = -c1%Im(:,nh)*ad_k1a1(k)
         z%Im(:,nh) =  c1%Re(:,nh)*ad_k1a1(k)
      _loop_km_end

   end subroutine var_coll_grad


!------------------------------------------------------------------------
!  find the divergence of a vector
!------------------------------------------------------------------------
   subroutine var_coll_div(r,t,z, dv)
      type (coll), intent(in)  :: r,t,z
      type (coll), intent(out) :: dv
      _loop_km_vars

      call var_coll_meshmult(1,mes_D%dr(1),r, c1)
      _loop_km_begin
         c2%Re(:,nh) =  &
	    r%Re(:,nh)*mes_D%r(:,-1) + c1%Re(:,nh)  &
	  - t%Im(:,nh)*ad_m1r1(:,m) - z%Im(:,nh)*ad_k1a1(k)
         c2%Im(:,nh) =  &
	    r%Im(:,nh)*mes_D%r(:,-1) + c1%Im(:,nh)  &
	  + t%Re(:,nh)*ad_m1r1(:,m) + z%Re(:,nh)*ad_k1a1(k)
      _loop_km_end
      call var_coll_copy(c2, dv)
      
   end subroutine var_coll_div


!------------------------------------------------------------------------
!  shift solution in z and theta
!------------------------------------------------------------------------
   subroutine var_coll_shift(dt,dz,a, b)
      double precision, intent(in)  :: dt,dz
      type (coll),      intent(in)  :: a
      type (coll),      intent(out) :: b
      double precision, save :: sink(-i_K1:i_K1), sinm(0:i_M1), dz_=1d8
      double precision, save :: cosk(-i_K1:i_K1), cosm(0:i_M1), dt_=1d8
      double precision :: dRe(i_N)
      _loop_km_vars
      
      call var_coll_copy(a, b)

      if(dz/=0d0) then
         if(dz/=dz_) then
            dz_ = dz      
            do k = -i_K1, i_K1
               sink(k) = dsin(dz*k*d_alpha)
               cosk(k) = dcos(dz*k*d_alpha)
            end do 
         end if
         _loop_km_begin
            dRe = b%Re(:,nh)
            b%Re(:,nh) = dRe*cosk(k) + b%Im(:,nh)*sink(k)
            b%Im(:,nh) = b%Im(:,nh)*cosk(k) - dRe*sink(k)
         _loop_km_end
      end if
      
      if(dt/=0d0) then
         if(dt/=dt_) then
            dt_ = dt      
            do m = 0, i_M1
               sinm(m) = dsin(dt*m*i_Mp)
               cosm(m) = dcos(dt*m*i_Mp)
            end do 
         end if
         _loop_km_begin
            dRe = b%Re(:,nh)
            b%Re(:,nh) = dRe*cosm(m) + b%Im(:,nh)*sinm(m)
            b%Im(:,nh) = b%Im(:,nh)*cosm(m) - dRe*sinm(m)
         _loop_km_end
      end if
   
   end subroutine var_coll_shift


!------------------------------------------------------------------------
!   i==0:  divergence-free conditions only
!   i==1:  mirror symmetric (Z)
!   i==2:  shift & reflect (S)
!   i==3:  shift & rotate (Omega)
!   NB: 'highly symm' == sh&ref + sh&rot
!------------------------------------------------------------------------
 subroutine var_imposesym(i, r,t,z)   
   integer,    intent(in)    :: i
   type(coll), intent(inout) :: r,t,z
   double precision :: a, b, bsend(i_N,6)
   integer :: nh_,p1,p2
   _loop_km_vars

   if(i==1 .or. i==2) then
      do m = 0, i_M1
         do k = -i_K1, i_K1
            if(k<0 .and. m==0) cycle
            
            call var_harmonic( k,m, nh ,p1)
            call var_harmonic(-k,m, nh_,p2)
            if(m==0) nh_=nh
            if(m==0) p2 =p1

            if(mpi_rnk==p1) then
               b = -1d0
               a =  1d0
               if(i==2 .and. modulo(k,2)/=0) a = -1d0
               if(m==0) b = 1d0
               bsend(:,1) =  r%Re(:,nh) * a
               bsend(:,2) =  r%Im(:,nh) * a * b
               bsend(:,3) = -t%Re(:,nh) * a
               bsend(:,4) = -t%Im(:,nh) * a * b
               bsend(:,5) =  z%Re(:,nh) * a
               bsend(:,6) =  z%Im(:,nh) * a * b
            end if
#ifdef _MPI
            if(p1/=p2) then
               mpi_tg = m*(2*i_K1+1) + k
               if(mpi_rnk==p1)  &
                  call mpi_send( bsend, i_N*6, mpi_double_precision,  &
                     p2, mpi_tg, mpi_comm_world, mpi_er)
               if(mpi_rnk==p2)  &
                  call mpi_recv( bsend, i_N*6, mpi_double_precision,  &
                     p1, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            end if
#endif         
            if(mpi_rnk==p2) then
               c1%Re(:,nh_) = bsend(:,1)
               c1%Im(:,nh_) = bsend(:,2)
               c2%Re(:,nh_) = bsend(:,3)
               c2%Im(:,nh_) = bsend(:,4)
               c3%Re(:,nh_) = bsend(:,5)
               c3%Im(:,nh_) = bsend(:,6)
            end if
         end do
      end do

      _loop_km_begin
         r%Re(:,nh) = ( r%Re(:,nh) + c1%Re(:,nh) ) / 2d0
         r%Im(:,nh) = ( r%Im(:,nh) + c1%Im(:,nh) ) / 2d0
         t%Re(:,nh) = ( t%Re(:,nh) + c2%Re(:,nh) ) / 2d0
         t%Im(:,nh) = ( t%Im(:,nh) + c2%Im(:,nh) ) / 2d0
         z%Re(:,nh) = ( z%Re(:,nh) + c3%Re(:,nh) ) / 2d0
         z%Im(:,nh) = ( z%Im(:,nh) + c3%Im(:,nh) ) / 2d0
      _loop_km_end
   end if
   
   if(i==3) then
      _loop_km_begin
         if(modulo(k+m,2)==0) cycle
         r%Re(:,nh) = 0d0
         r%Im(:,nh) = 0d0
         t%Re(:,nh) = 0d0
         t%Im(:,nh) = 0d0
         z%Re(:,nh) = 0d0
         z%Im(:,nh) = 0d0
      _loop_km_end
   end if
   
   if(mpi_rnk==0) then
      r%Re(:,0) = 0d0
      r%Im(:,0) = 0d0
      t%Im(:,0) = 0d0
      z%Im(:,0) = 0d0
   end if

 end subroutine var_imposesym 


!------------------------------------------------------------------------
!  get index of nearest point
!------------------------------------------------------------------------
   subroutine var_nearestpt(x, ix)
      double precision, intent(in)  :: x(3)
      integer,          intent(out) :: ix(3)
      double precision, save :: dt,dz, r_(i_N)
      logical, save :: set = .false.

      if(.not.set) then 
         r_ = mes_D%r(:,1)
         dt = 2d0*d_PI/dble(i_Mp*i_Th)
         dz = 2d0*d_PI/(d_alpha*i_Z)
         set = .true.
      end if
      if(x(1)<0d0 .or. x(1)>1d0)              stop 'var_nearestpt: err1'
      if(x(2)<0d0 .or. x(2)>2d0*d_PI)         stop 'var_nearestpt: err2'
      if(x(3)<0d0 .or. x(3)>2d0*d_PI/d_alpha) stop 'var_nearestpt: err3'

      ix(1) = 1
      do while( dabs(x(1)-r_(ix(1))) > dabs(x(1)-r_(ix(1)+1)) )
         ix(1) = ix(1) + 1
         if(ix(1)==i_N) exit
      end do
      ix(2) = modulo(nint(x(2)/dt),i_Th)
      ix(3) = modulo(nint(x(3)/dz),i_Z)
      
   end subroutine var_nearestpt


!------------------------------------------------------------------------
!  trilinear interpolation about the nearest point
!------------------------------------------------------------------------
   subroutine var_interpnearpt(x,ix,p, c)
      double precision, intent(in)  :: x(3)
      integer,          intent(in)  :: ix(3)
      type (phys),      intent(in)  :: p
      double precision, intent(out) :: c
      double precision, save :: dt,dz, r_(i_N)
      logical, save :: set = .false.
      double precision :: rd,td,zd, c00,c10,c01,c11, c0,c1
      integer :: ir0,ir1,it0,it1,iz0,iz1

      if(.not.set) then 
         r_ = mes_D%r(:,1)
         dt = 2d0*d_PI/dble(i_Mp*i_Th)
         dz = 2d0*d_PI/(d_alpha*i_Z)
         set = .true.
      end if

      if(ix(1)==1 .or. x(1)>r_(ix(1))) then
         ir0 = ix(1)
         ir1 = ix(1)+1
      else
         ir0 = ix(1)-1
         ir1 = ix(1)
      end if
      rd = (x(1)-r_(ir0))/(r_(ir1)-r_(ir0))

      it0 = nint(x(2)/dt)/i_Th
      td = x(2) - dt*(ix(2) + i_Th*it0)
      if(td<0d0) then
         td  = td/dt + 1d0
         it0 = modulo(ix(2)+i_Th-1,i_Th)
         it1 = ix(2)
      else
         td  = td/dt
         it0 = ix(2)
         it1 = modulo(ix(2)+1,i_Th)
      end if

      iz0 = nint(x(3)/dz)/i_Z
      zd = x(3) - dz*(ix(3) + i_Z*iz0)
      if(zd<0d0) then
         zd  = zd/dz + 1d0
         iz0 = modulo(ix(3)+i_Z-1,i_Z)
         iz1 = ix(3)
      else
         zd  = zd/dz
         iz0 = ix(3)
         iz1 = modulo(ix(3)+1,i_Z)
      end if

      c00 = p%Re(iz0,it0,ir0)*(1d0-rd) + p%Re(iz0,it0,ir1)*rd
      c10 = p%Re(iz0,it1,ir0)*(1d0-rd) + p%Re(iz0,it1,ir1)*rd
      c01 = p%Re(iz1,it0,ir0)*(1d0-rd) + p%Re(iz1,it0,ir1)*rd
      c11 = p%Re(iz1,it1,ir0)*(1d0-rd) + p%Re(iz1,it1,ir1)*rd
      c0 = c00*(1d0-td) + c10*td
      c1 = c01*(1d0-td) + c11*td
      c = c0*(1d0-zd) + c1*zd      
      
 end subroutine var_interpnearpt


!------------------------------------------------------------------------
!  norm   (1/2) \int a.a dV  for each l,m
!------------------------------------------------------------------------
   subroutine var_coll_norm(a, E,Ek,Em)
      type (coll),      intent(in)  :: a
      double precision, intent(out) :: E, Ek(0:i_K1), Em(0:i_M1)
#ifdef _MPI
      double precision :: E_, Ek_(0:i_K1), Em_(0:i_M1)
#endif
      double precision :: w, b, f(i_N)
      _loop_km_vars
      
      Ek = 0d0
      Em = 0d0
      w  = 4d0 * d_PI*d_PI / d_alpha
      if(d_alpha==0d0) stop 'var_coll_norm: alpha=0'
      _loop_km_begin
         f =  w * ( a%Re(:,nh)*a%Re(:,nh)  &
                  + a%Im(:,nh)*a%Im(:,nh) )
         b = dot_product(f,mes_D%intrdr )
         if(k==0 .and. m==0) b = 0.5d0*b 
         Em(m)      = Em(m)      + b
         Ek(abs(k)) = Ek(abs(k)) + b
      _loop_km_end

#ifdef _MPI
      call mpi_allreduce( Ek, Ek_, 1+i_K1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Ek = Ek_
      call mpi_allreduce( Em, Em_, 1+i_M1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em = Em_
#endif
      
      E = sum(Em)
   
   end subroutine var_coll_norm

!*************************************************************************
 end module variables
!*************************************************************************
