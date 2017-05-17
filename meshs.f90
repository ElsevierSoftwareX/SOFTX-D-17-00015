!**************************************************************************
! Definition of mesh-types.
!
! Public variables:
!
!   r(n,k)              radial points r_j^k, j in [1,N]
!
!   dr(p)%M             p-th derivative of f  =  dr(p)%M . f
!
!   radLap%M            radial part of laplacian  (1/r)dr + drr
!
!
! (mesh)-type matricies are stored s.t.   M(   KL+1+n-j, j) = A(n,j)
! (lumesh)-type matricies are stored s.t. M( 2*KL+1+n-j, j) = A(n,j)
!
! Note that the index n is usually used to denote r_n and the operator on
! the n-th point corresponds to the n-th row of the matrix.  The index
! j usually denotes the column and is therefore in the outer loop below.
! Matrix multiplication  y := M . x  can be achieved accessing M only
! once along the memory by the loop form
!
!   y = 0
!   do j = 1, N
!      do n = max(1,j-KL), min(j+KL,N)
!         y(n) = y(n) + M(KL+1+n-j,j) * x(j)
!
!
!**************************************************************************
#include "../parallel.h"
 module meshs
!**************************************************************************
   use mpif
   use parameters
   implicit none
   save
                                ! M(KL+1+n-j, j) = A(n,j)
   type mesh
      double precision :: M(2*i_KL+1, i_N)
   end type mesh
                                ! M(2*KL+1+n-j, j) = A(n,j)
                                ! see lapack dgbtrf
   type lumesh
      integer          :: ipiv(i_N)
      double precision :: M(3*i_KL+1, i_N)
   end type lumesh


   type rdom
      integer          :: N, pNi,pN, pNi_(0:_Np-1),pN_(0:_Np-1)
      double precision :: r(i_N,-3:3)
      double precision :: intrdr(i_N)
      double precision :: dr0(1:1+i_KL,0:i_KL), dr1(1:1+i_KL,0:i_KL)
      type (mesh)      :: dr(i_KL)
      type (mesh)      :: radLap
   end type rdom


   type (rdom) :: mes_D

 contains


!------------------------------------------------------------------------
!  Precompute the mesh stuff
!------------------------------------------------------------------------
   subroutine mes_precompute()
      double precision :: dr
      integer :: n, N_
      logical :: file_exist
         			! collocation pts
                                ! T_{N-1}(x) has N extrema
                                ! Numerical Recipes 5.8.5
      if(_Nr>i_N) stop '_Nr>i_N'
      mes_D%N = i_N

      inquire(file='mesh.in', exist=file_exist)
      if(file_exist) then
         if(mpi_rnk==0) print*, ' using custom mesh!'
         open(99,status='old',file='mesh.in')
         read(99,*) n
         if(n/=i_N) stop 'mes_precompute: mesh.in err1: N mismatch'
         do n = 1, i_N
            read(99,*) mes_D%r(n,1) 
         end do
         close(99)
         if(mes_D%r(1,1)<0d0) stop 'mes_precompute: mesh.in err 2'
         dr = dabs(mes_D%r(i_N,1)-1d0)
         if(dr>1d-8) stop 'mes_precompute: mesh.in err 3'
         
      else ! default mesh
         N_ = i_N+int(dsqrt(dble(i_N)))
         do n = N_-i_N+1, N_
            mes_D%r(n-N_+i_N,1) = 0.5d0*( 1d0+dcos(d_PI*(N_-n)/dble(N_)) )
         end do
         do n = 1, 10
            dr = 1.5d0*mes_D%r(1,1)-0.5d0*mes_D%r(2,1)
            mes_D%r(:,1) = mes_D%r(:,1)*(1d0+dr) - dr
         end do
      end if
      
      call mes_rdom_init(mes_D)

   end subroutine mes_precompute


!------------------------------------------------------------------------
! for the radial domain D, set the number of points D%N and their
! values D%r, then mes_rdom_init() will fill in the rest...
!------------------------------------------------------------------------
   subroutine mes_rdom_init(D)
      type (rdom), intent(inout) :: D
      integer, parameter :: lda = 2*i_KL+1
      double precision   :: A(lda,lda), w(lda), c,e
      double precision   :: r_(0:D%N)
      integer :: n,i,j,nn, l,r
         			! distribute radial points over procs
      D%pN_ = 0
      do n = 1, D%N
         r = _Nr - modulo(n-1,_Nr) - 1
         D%pN_(r) = D%pN_(r) + 1
      end do
      D%pNi_(0) = 1
      do r = 1, _Nr-1
         D%pNi_(r) = D%pNi_(r-1) + D%pN_(r-1)
      end do
      do r = _Nr, _Np-1
         D%pN_ (r) = D%pN_ (r-_Nr)
         D%pNi_(r) = D%pNi_(r-_Nr)
      end do
      D%pN  = D%pN_ (mpi_rnk)
      D%pNi = D%pNi_(mpi_rnk)
         			! get ad_r(:,j) = r^j
      do j = -3, 3
         if(j==1) cycle
         D%r(:,j) = D%r(:,1)**j
      end do

                                ! get finite difference weights    
      if(i_KL<2) stop 'mes: KL<2'
      do i = 1, i_KL
         D%dr(i)%M = 0d0
      end do
      D%radLap%M = 0d0
      D%intrdr = 0d0

      r_(1:D%N) = D%r(:,1)
      r_(0) = 0d0

      do n = 1, D%N
         l = max(1,n-i_KL)
         r = min(n+i_KL,D%N)
         nn = r-l+1
	 do i = 1, i_KL
   	    call mes_weights(i,nn,r_(l),r_(n), w)
            do j = l, r
               D%dr(i)%M(i_KL+1+n-j,j) = w(j-l+1)
            end do
	 end do
      end do
         			! weights for extrapolation to r=0, 1
      do i = 0, i_KL
         call mes_weights(i,1+i_KL,r_(1),r_(0), D%dr0(1,i))
         call mes_weights(i,1+i_KL,r_(D%N-i_KL),r_(D%N), D%dr1(1,i))
      end do
                                ! get  (1/r)dr + drr
      do j = 1, D%N
         l = max(1,j-i_KL)
         r = min(j+i_KL,D%N)
         do n = l, r
            D%radLap%M(i_KL+1+n-j,j)  &
               = D%dr(2)%M(i_KL+1+n-j,j)  &
               + D%dr(1)%M(i_KL+1+n-j,j) / r_(n)
         end do
      end do
            			! weights for integration  r dr
      do n = 1, D%N
         l = max(1,n-i_KL)
         r = min(n+i_KL,D%N)
         nn = r-l+1
	 do i = 1, nn
	    call mes_weights(i-1,nn,r_(l),r_(n), w)
	    A(i,:) = w 
	 end do
         c = 1d0
         e = 1d0
         do i = 1, nn
            c = c * (r_(min(n+1,D%N))-r_(n)) / dble(i)
            e = e * (r_(max(  0,n-1))-r_(n)) / dble(i)
            D%intrdr(l:l+nn-1) =  &
               D%intrdr(l:l+nn-1) + 0.5d0*(c-e)*A(i,1:nn)*r_(n)
         end do
      end do

   end subroutine mes_rdom_init


!------------------------------------------------------------------------
!  given n pts x, find weights w to get i^th deriv at x0
!------------------------------------------------------------------------
   subroutine mes_weights(i,n,x,x0, w)
      integer,          intent(in)  :: i, n
      double precision, intent(in)  :: x(n), x0
      double precision, intent(out) :: w(n)
      double precision :: A(n,n), d(n), x_(n), x0_, xl, xr
      integer :: j
   
      xl = min(x0,minval(x))
      xr = max(x0,maxval(x))
      x_  = (x  - xl)/(xr-xl)
      x0_ = (x0 - xl)/(xr-xl)
      do j = 1, n
         call cheby(j-1, 0,  x_, n, A(1,j))
         call cheby(j-1, i, x0_, 1, d(j))
      end do
      call mes_mat_invert(n,A,n)
      w = matmul(d,A) / ((xr-xl)**i)
/*
      A(:,1) = 1d0
      do j = 2, n
         A(:,j) = A(:,j-1) * (x-x0) / dble(j-1)
      end do
      call mes_mat_invert(n,A,n)
      w = A(i+1,:)
*/	 
   end subroutine mes_weights 
	 
	 
!------------------------------------------------------------------------
!  Replace nxn matrix A by its inverse
!------------------------------------------------------------------------
   subroutine mes_mat_invert(n,A,lda)
      integer, intent(in) :: n, lda
      double precision, intent(inout) :: A(lda,n)
      double precision :: work(n)
      integer :: info, ipiv(n)
   
      call dgetrf(n, n, A, lda, ipiv, info)
      if(info /= 0) stop 'matrix inversion error 1'

      call dgetri(n, A, lda, ipiv, work, n, info)
      if(info /= 0) stop 'matrix inversion error 2'

   end subroutine mes_mat_invert

!------------------------------------------------------------------------
!  Replace complex nxn matrix A by its inverse
!------------------------------------------------------------------------
 subroutine mes_mat_invert_complex(n,Ar,Ai,lda)
      integer, intent(in) :: n, lda
      double precision, intent(inout) :: Ar(lda,n), Ai(lda,n)
      double complex :: A(lda,n), work(n)
      integer :: info, ipiv(n)
   
      A = dcmplx(Ar,Ai)

      call zgetrf(n, n, A, lda, ipiv, info)
      if(info /= 0) stop 'mat_invert_complex, zgetrf fail'

      call zgetri(n, A, lda, ipiv, work, n, info)
      if(info /= 0) stop 'mat_invert_complex, zgetri fail'
      
      Ar =  dble(A)
      Ai = dimag(A)

 end subroutine mes_mat_invert_complex

!**************************************************************************
 end module meshs
!**************************************************************************
