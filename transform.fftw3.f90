!*************************************************************************
! Transformation using the FFT, de-aliased.
! Written for FFTW version 3 -- documentation at  www.fftw.org
! See FFTW reference, sections 2.3, 4.3.3, 4.4.1-2
!
! \sum_{-Ma}^{Ma} == \sum_0^{3M-1}
! logical size of DFT n==3M=Th
!
!*************************************************************************
#include "../parallel.h"
 module transform
!*************************************************************************
   use parameters
   use variables
   implicit none
   save

   integer, parameter, private :: i_3K = 3*i_K
   integer, parameter, private :: i_3M = 3*i_M
   integer, parameter, private :: i_Ma = (3*i_M)/2
   double complex,     private :: X(0:i_3K-1, 0:_Ms-1)
   double complex,     private :: Y(0:i_3K-1, 0:_Ms-1)
   double complex,     private :: Xs(0:i_pZ-1, 0:i_Ma)
   double precision,   private :: Ys(0:i_pZ-1, 0:i_3M-1)
   integer*8,          private :: plan_c2cf, plan_c2cb, plan_r2c, plan_c2r

   double complex,     private :: T(0:i_3K-1, 0:_Ms-1, i_pN)
   double complex,     private :: Ts(0:i_pZ-1, 0:i_M1, i_pN)

 contains

!------------------------------------------------------------------------
! Setup plans for transforms.  
!------------------------------------------------------------------------
   subroutine tra_precompute()
      integer, parameter :: flag=32 !=FFTW_PATIENT see fftw3.f
      integer :: sgn, n(1), howmany, inembed(1), onembed(1)
 
      n = (/i_3K/)            
      howmany = _Ms
      inembed = (/i_3K*_Ms/)
      onembed = (/i_3K*_Ms/)
      sgn = 1
      call dfftw_plan_many_dft(plan_c2cf, 1, n, howmany,  &
         X, inembed, 1, i_3K,  Y, onembed, 1, i_3K,  sgn, flag)
      sgn = -1
      call dfftw_plan_many_dft(plan_c2cb, 1, n, howmany,  &
         Y, onembed, 1, i_3K,  X, inembed, 1, i_3K,  sgn, flag)

      n = (/i_3M/)
      howmany = i_pZ
      inembed = (/i_pZ*(i_Ma+1)/)
      onembed = (/i_pZ*i_3M/)
      call dfftw_plan_many_dft_c2r(plan_c2r, 1, n, howmany,  &
         Xs, inembed, i_pZ, 1,  Ys, onembed, i_pZ, 1,  flag)
      call dfftw_plan_many_dft_r2c(plan_r2c, 1, n, howmany,  &
         Ys, onembed, i_pZ, 1,  Xs, inembed, i_pZ, 1,  flag)
      
   end subroutine tra_precompute


!------------------------------------------------------------------------
!  Convert spectral to real space
!------------------------------------------------------------------------
   subroutine tra_spec2phys(s, p)
      type (spec), intent(in)  :: s
      type (phys), intent(out) :: p
      integer :: nh, n,m,m_
      				! for each r_n ...      
      do n = 1, mes_D%pN
         if(mpi_rnk/_Nr==0) then
            X(0:i_K1,   0) = dcmplx(s%Re(0:i_K1,n),s%Im(0:i_K1,n))
            X(i_K:2*i_K,0) = 0d0
            X(2*i_K+1:, 0) = dcmplx(s%Re(i_K1:1:-1,n),-s%Im(i_K1:1:-1,n))
            m_ = 1
            nh = 2*i_K-1
         else
            m_ = 0
            nh = i_K1
         end if
         do m = m_, _Ms1
            X(0:i_K1,   m) = dcmplx(s%Re(nh:nh+i_K1,n),s%Im(nh:nh+i_K1,n))
            X(i_K:2*i_K,m) = 0d0
            X(2*i_K+1:, m) = dcmplx(s%Re(nh-i_K1:nh-1,n),s%Im(nh-i_K1:nh-1,n))
            nh = nh + 2*i_K-1
         end do
         call dfftw_execute(plan_c2cf)
#if _Ns == 1 
         Xs(:,0:i_M1) = Y
#else
         T(:,:,n) = Y
      end do
      call tra_T2Ts()
      do n = 1, mes_D%pN
         Xs(:,0:i_M1) = Ts(:,:,n)         
#endif
         Xs(:,i_M:) = 0d0
         call dfftw_execute(plan_c2r)
         p%Re(:,:,n) = Ys
      end do
   
   end subroutine tra_spec2phys


!------------------------------------------------------------------------
!  Convert real to spectral space
!------------------------------------------------------------------------
   subroutine tra_phys2spec(p, s)
      type (phys), intent (in)  :: p
      type (spec), intent (out) :: s
      integer :: nh, n,m, m_
      double precision :: scale_
         			! scale, FFTW 4.7.2
      scale_ = 1d0 / dble(i_3K*i_3M)

      do n = 1, mes_D%pN
         Ys = scale_ * p%Re(:,:,n)
         call dfftw_execute(plan_r2c)
#if _Ns == 1
         Y = Xs(:,0:i_M1)
#else
         Ts(:,:,n) = Xs(:,0:i_M1)
      end do
      call tra_Ts2T()
      do n = 1, mes_D%pN
         Y = T(:,:,n)
#endif
         call dfftw_execute(plan_c2cb)
         if(mpi_rnk/_Nr==0) then
            s%Re(0:i_K1,n) =  dble(X(0:i_K1,0))
            s%Im(0:i_K1,n) = dimag(X(0:i_K1,0))
            m_ = 1
            nh = 2*i_K-1
         else
            m_ = 0
            nh = i_K1
         end if
         do m = m_, _Ms1
            s%Re(nh:nh+i_K1,n) =  dble(X(0:i_K1,m))
            s%Im(nh:nh+i_K1,n) = dimag(X(0:i_K1,m))
            s%Re(nh-i_K1:nh-1,n) =  dble(X(2*i_K+1:,m))
            s%Im(nh-i_K1:nh-1,n) = dimag(X(2*i_K+1:,m))
            nh = nh + 2*i_K-1
         end do
      end do
         
   end subroutine tra_phys2spec


!------------------------------------------------------------------------
! transposes
!------------------------------------------------------------------------
#if _Ns != 1
   subroutine tra_T2Ts()
      double precision :: bsend(2*i_pN*_Ms*i_pZ,0:_Ns-1)
      double precision :: brecv(2*i_pN*_Ms*i_pZ,0:_Ns-1)
      integer :: stp, dst,src, l,j, rnk,rko
      integer :: n,m, pm0,jz0

      rnk = mpi_rnk/_Nr
      rko = modulo(mpi_rnk,_Nr)

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*mes_D%pN*_Ms*i_pZ,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         dst  = modulo(stp+rnk, _Ns)*_Nr + rko
         jz0  = (dst/_Nr)*i_pZ
         l = 1
         do n = 1, mes_D%pN
            do m = 0, _Ms1
               do j = jz0, jz0+i_pZ-1
                  bsend(l,  stp) =  dble(T(j,m,n))
                  bsend(l+1,stp) = dimag(T(j,m,n))
                  l = l + 2
               end do
            end do
         end do
         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*mes_D%pN*_Ms*i_pZ,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do
 
      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         pm0 = (src/_Nr)*_Ms
         l = 1
         do n = 1, mes_D%pN
            do m = pm0, pm0+_Ms1
               do j = 0, i_pZ-1
                  Ts(j,m,n) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
            end do
         end do
      end do

      do stp = 0, _Ns-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do

   end subroutine tra_T2Ts

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   subroutine tra_Ts2T()
      double precision :: bsend(2*i_pN*_Ms*i_pZ,0:_Ns-1)
      double precision :: brecv(2*i_pN*_Ms*i_pZ,0:_Ns-1)
      integer :: stp, dst,src, l,j, rnk,rko
      integer :: n,m, pm0,jz0

      rnk = mpi_rnk/_Nr
      rko = modulo(mpi_rnk,_Nr)

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         mpi_tg = stp + rko
         call mpi_irecv( brecv(1,stp), 2*i_pN*_Ms*i_pZ,  &
            mpi_double_precision, src, mpi_tg, mpi_comm_world,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         dst  = modulo(stp+rnk, _Ns)*_Nr + rko
         pm0 = (dst/_Nr)*_Ms
         l = 1
         do n = 1, mes_D%pN
            do m = pm0, pm0+_Ms1
               do j = 0, i_pZ-1
                  bsend(l,  stp) =  dble(Ts(j,m,n))
                  bsend(l+1,stp) = dimag(Ts(j,m,n))
                  l = l + 2
               end do
            end do
         end do
         mpi_tg = stp + rko
         call mpi_isend( bsend(1,stp), 2*i_pN*_Ms*i_pZ,  &
            mpi_double_precision, dst, mpi_tg, mpi_comm_world,  &
            mpi_rq(mpi_sze+stp), mpi_er)
      end do

      do stp = 0, _Ns-1
         src  = modulo(mpi_sze-stp+rnk, _Ns)*_Nr + rko
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         jz0  = (src/_Nr)*i_pZ
         l = 1
         do n = 1, mes_D%pN
            do m = 0, _Ms1
               do j = jz0, jz0+i_pZ-1
                  T(j,m,n) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
            end do
         end do
      end do

      do stp = 0, _Ns-1
         call mpi_wait( mpi_rq(mpi_sze+stp), mpi_st, mpi_er)
      end do

   end subroutine tra_Ts2T
#endif

!*************************************************************************
 end module transform
!*************************************************************************
