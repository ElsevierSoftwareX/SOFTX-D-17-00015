!*************************************************************************
#include "../parallel.h"
 module dx
!*************************************************************************
   use io
   implicit none
   save

   integer :: var,ddom,dvec,zshft,submean
   double precision :: zval,zval2,thval,rval, zshfta,zshftb
   double precision :: zoffset,toffset,velprof(i_N)

   double precision :: dtprev, d_Re_
   double precision :: maxA(3),minA(3),maxA_(3),minA_(3)
   character(4) :: cnum
   integer :: fl,fls,ffrst,flast,fstep
   type (phys), private :: pr,pt,pz
   type (coll), private :: c1,c2,c3,c4,c5,c6
   type (spec), private :: s1,s2,s3
 
   double precision, private :: GT(3,3), ST(3,3), AT(3,3), E_(3)
   type (phys), private :: p21,p22,p23,p31,p32,p33

   integer :: mat_Ads, mat_nx, mat_ny, mat_nz
   double precision, allocatable, private :: mat_A(:,:,:,:)
   double precision, allocatable, private :: mat_x(:) 
   double precision, allocatable, private :: mat_y(:) 
   double precision, allocatable, private :: mat_z(:)

 contains 
 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_precompute()
   end subroutine dx_precompute


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_preparedata()
      integer :: e,f,n 
      double precision :: d, a(i_N), b(i_N)
      _loop_km_vars

      vel_uz%Re(:,0) = vel_uz%Re(:,0) + vel_U
      if(submean==1) then
         vel_uz%Re(:,0) = vel_uz%Re(:,0) - vel_U
      else if(submean==2) then
         vel_uz%Re(:,0) = 0d0
      else if(submean==3) then
         vel_uz%Re(:,0) = vel_uz%Re(:,0) - velprof
      end if

      if(zshft/=0 .or. toffset/=0d0) then
         d = 0d0
         if(zshft==1) then
            stop 'puf_ module not implemented'
         else if(zshft==2) then
            if(d_time>=0d0) stop 'set d_time<0d0'
            d = -(zshfta*tim_t + zshftb) + zoffset
         else if(zshft==3) then
            n = -1
            do while(n/=fl)
               read(77,*) n, d
            end do
            d = -d + zoffset
         end if
         call var_coll_shift(toffset,d,vel_ur, vel_ur)
         call var_coll_shift(toffset,d,vel_ut, vel_ut)
         call var_coll_shift(toffset,d,vel_uz, vel_uz)
      end if


      if(var==1) then
         call var_coll2spec(vel_ur,s1, c2=vel_ut,s2=s2, c3=vel_uz,s3=s3)
         call tra_spec2phys( s1, pr)
         call tra_spec2phys( s2, pt)
         call tra_spec2phys( s3, pz)
      
      else if(var==2) then
         call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
         call var_coll2spec(c1,s1, c2=c2,s2=s2, c3=c3,s3=s3)
         call tra_spec2phys( s1, pr)
         call tra_spec2phys( s2, pt)
         call tra_spec2phys( s3, pz)

      else if(var==3) then
         call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
         call var_coll2spec(c3, s1)
         call tra_spec2phys(s1, pr)

      else if(var==4) then
         call var_coll_curl(vel_ur,vel_ut,vel_uz, c1,c2,c3)
         call var_coll_curl(c1,c2,c3, c1,c2,c3)
         call var_coll2spec(c3, s1)
         call tra_spec2phys(s1, pr)

      else if(var==5) then
         vel_uz%Re(:,1:) = 0d0
         vel_uz%Im(:,1:) = 0d0
         call var_coll2spec(vel_uz, s1)
         call tra_spec2phys(s1, pr)

      else if(var==6) then
         stop ' option 6 not implemented'
!         vel_P2%Re(:,0) = 0d0
!         vel_P2%Im(:,0) = 0d0
!         call var_coll2spec(vel_P2, sp)
!         call tra_spec2phys(sp, pr)
      
      else if(var==7) then
         stop ' option 7 not implemented'
         
      else if(var==8) then

         call var_coll_copy(vel_ur, c1)
         call var_coll_copy(vel_ut, c2)
         call var_coll_copy(vel_uz, c3)
         if(submean/=0)  &
            print*, 'Warning: only strictly valid for total velocity.'
                                      
         call var_coll_meshmult(1,mes_D%dr(1),c1,c4)
         call var_coll_meshmult(1,mes_D%dr(1),c2,c5)
         call var_coll_meshmult(0,mes_D%dr(1),c3,c6)

         call var_coll2spec(c4,s1, c2=c5,s2=s2, c3=c6,s3=s3)
         call tra_spec2phys( s1, pr)
         call tra_spec2phys( s2, pt)
         call tra_spec2phys( s3, pz)
         
         _loop_km_begin
            c4%Re(:,nh) = mes_D%r(:,-1)*(-c1%Im(:,nh)*m*i_Mp-c2%Re(:,nh))
            c4%Im(:,nh) = mes_D%r(:,-1)*( c1%Re(:,nh)*m*i_Mp-c2%Im(:,nh))
            c5%Re(:,nh) = mes_D%r(:,-1)*(-c2%Im(:,nh)*m*i_Mp+c1%Re(:,nh))
            c5%Im(:,nh) = mes_D%r(:,-1)*( c2%Re(:,nh)*m*i_Mp+c1%Im(:,nh))
            c6%Re(:,nh) = mes_D%r(:,-1)*(-c3%Im(:,nh)*m*i_Mp)
            c6%Im(:,nh) = mes_D%r(:,-1)*( c3%Re(:,nh)*m*i_Mp)
         _loop_km_end
         call var_coll2spec(c4,s1, c2=c5,s2=s2, c3=c6,s3=s3)
         call tra_spec2phys( s1, p21)
         call tra_spec2phys( s2, p22)
         call tra_spec2phys( s3, p23)
      
         _loop_km_begin
            c4%Re(:,nh) = -c1%Im(:,nh)*d_alpha*k
            c4%Im(:,nh) =  c1%Re(:,nh)*d_alpha*k
            c5%Re(:,nh) = -c2%Im(:,nh)*d_alpha*k
            c5%Im(:,nh) =  c2%Re(:,nh)*d_alpha*k
            c6%Re(:,nh) = -c3%Im(:,nh)*d_alpha*k
            c6%Im(:,nh) =  c3%Re(:,nh)*d_alpha*k
         _loop_km_end
         call var_coll2spec(c4,s1, c2=c5,s2=s2, c3=c6,s3=s3)
         call tra_spec2phys( s1, p31)
         call tra_spec2phys( s2, p32)
         call tra_spec2phys( s3, p33)

         do n = 1, i_N
            do m = 0, i_Th-1
               do k = 0, i_Z-1
                  GT(1,1) = pr%Re(k,m,n)
                  GT(1,2) = pt%Re(k,m,n)
                  GT(1,3) = pz%Re(k,m,n)
                  GT(2,1) = p21%Re(k,m,n)
                  GT(2,2) = p22%Re(k,m,n)
                  GT(2,3) = p23%Re(k,m,n)
                  GT(3,1) = p31%Re(k,m,n)
                  GT(3,2) = p32%Re(k,m,n)
                  GT(3,3) = p33%Re(k,m,n)
                  ST = 0.5d0 * (GT + transpose(GT))
                  AT = 0.5d0 * (GT - transpose(GT))
                  GT = matmul(ST,ST) + matmul(AT,AT)
                  call get_eigvals(3,GT,3,E_)             
                  pr%Re(k,m,n) = E_(2)               
               end do
            end do
         end do      

      end if   

            
   end subroutine dx_preparedata


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_volume()
      double precision :: dx,dy,dz, r,t,z,a(3),x(3)
      integer :: n, i,j,k, ix(3)

      if(i_M==1) stop 'Need more theta points!  Use larger i_M.'
      if(i_K==1 .and. zval/=zval2) stop 'Need more z points! Use larger i_K.'

      allocate(mat_x(mat_nx))
      allocate(mat_y(mat_ny))
      allocate(mat_z(mat_nz))
      allocate(mat_A(mat_ny,mat_nx,mat_nz,mat_Ads))

      dx = 2d0/dble(mat_nx)
      dy = 2d0/dble(mat_ny)
      dz = (zval2-zval)/dble(mat_nz)
      
      do n = 1, mat_nx
         mat_x(n) = -1d0+dx*(-0.5d0+dble(n))
      end do
      do n = 1, mat_ny
         mat_y(n) = -1d0+dy*(-0.5d0+dble(n))
      end do
      do n = 1, mat_nz
         mat_z(n) = zval + dz*(n-1)
      end do

      do i = 1, mat_nx
         do j = 1, mat_ny
            do k = 1, mat_nz
               r = dsqrt(mat_x(i)**2+mat_y(j)**2)
               if(r==0d0) t = 0d0
               if(r/=0d0) t = atan2(mat_y(j),mat_x(i))   ! t in (-pi,pi]
               if(t<0d0)  t = t+2d0*d_PI   ! t in [0,2pi)
               z = mat_z(k)
               if(r>1d0) then
                  mat_A(j,i,k,:) = 0d0
               else
                  x = (/r,t,z/)
                  call var_nearestpt(x, ix)
                  if(mat_Ads==1) then
                     call var_interpnearpt(x,ix,pr, mat_A(j,i,k,1))
                  else  
                     call var_interpnearpt(x,ix,pr, a(1))
                     call var_interpnearpt(x,ix,pt, a(2))
                     call var_interpnearpt(x,ix,pz, a(3))
                     if(dvec==1) then
                        mat_A(j,i,k,1) = a(1)*cos(t) - a(2)*sin(t)
                        mat_A(j,i,k,2) = a(1)*sin(t) + a(2)*cos(t)
                        mat_A(j,i,k,3) = a(3)
                     else
                        mat_A(j,i,k,:) = a
                     end if
                  end if
               end if
            end do
         end do
      end do
     
   end subroutine dx_volume


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_rz_crsect()
      double precision :: dx, r,t,z,a(3),x(3)
      integer :: n,n_, i,j, ix(3)

      mat_nx = i_Z
      mat_ny = i_N*2
      mat_nz = 1
      allocate(mat_x(mat_nx))
      allocate(mat_y(mat_ny))
      allocate(mat_z(mat_nz))
      allocate(mat_A(mat_ny,mat_nx,mat_nz,mat_Ads))

      dx = (2d0*d_PI/d_alpha)/dble(mat_nx)
      do n = 1, mat_nx
         mat_x(n) = dx*(n-1)
      end do
      n_ = 0
      do n = -i_N, i_N
         if(n==0) cycle
         n_ = n_+1
         mat_y(n_) = mes_D%r(abs(n),1) * (n/abs(n))
      end do
      mat_z(1) = thval

      do i = 1, mat_nx
         do j = 1, mat_ny
            r = abs(mat_y(j))
            t = mat_z(1)
            if(mat_y(j)>0d0) t = t + d_PI
            if(t>2d0*d_PI) t = t - 2d0*d_PI
            z = mat_x(i)
            x = (/r,t,z/)
            if(i_M==1) x=(/r,0d0,z/)
            call var_nearestpt(x, ix)
            if(mat_Ads==1) then
               call var_interpnearpt(x,ix,pr, mat_A(j,i,1,1))
            else  
               call var_interpnearpt(x,ix,pr, a(1))
               call var_interpnearpt(x,ix,pt, a(2))
               call var_interpnearpt(x,ix,pz, a(3))
               if(dvec==1) then
                  mat_A(j,i,1,1) = a(1)*cos(t) - a(2)*sin(t)
                  mat_A(j,i,1,2) = a(1)*sin(t) + a(2)*cos(t)
                  mat_A(j,i,1,3) = a(3)
               else
                  mat_A(j,i,1,:) = a
               end if
            end if
         end do
      end do

   end subroutine dx_rz_crsect


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_thz_cylsect()
      double precision :: dx,dy, r,t,z,a(3),x(3)
      integer :: n, i,j, ix(3)

      mat_nx = i_Z
      mat_ny = i_Th*i_Mp
      mat_nz = 1
      allocate(mat_x(mat_nx))
      allocate(mat_y(mat_ny))
      allocate(mat_z(mat_nz))
      allocate(mat_A(mat_ny,mat_nx,mat_nz,mat_Ads))

      dx = (2d0*d_PI/d_alpha)/dble(mat_nx)
      do n = 1, mat_nx
         mat_x(n) = dx*(n-1)
      end do
      dy = 2d0*d_PI/dble(mat_ny)
      do n = 1, mat_ny
         mat_y(n) = dy*(n-1)
      end do
      mat_z(1) = rval

      do i = 1, mat_nx
         do j = 1, mat_ny
            r = rval
            t = mat_y(j)
            z = mat_x(i)
            x = (/r,t,z/)
            call var_nearestpt(x, ix)
            if(mat_Ads==1) then
               call var_interpnearpt(x,ix,pr, mat_A(j,i,1,1))
            else  
               call var_interpnearpt(x,ix,pr, a(1))
               call var_interpnearpt(x,ix,pt, a(2))
               call var_interpnearpt(x,ix,pz, a(3))
               if(dvec==1) then
                  mat_A(j,i,1,1) = a(1)*cos(t) - a(2)*sin(t)
                  mat_A(j,i,1,2) = a(1)*sin(t) + a(2)*cos(t)
                  mat_A(j,i,1,3) = a(3)
               else
                  mat_A(j,i,1,:) = a
               end if
            end if
         end do
      end do

   end subroutine dx_thz_cylsect


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_dump()
      integer :: e,f, ad,xd,yd,zd,dims(4),xid,yid,zid,Aid
      integer :: i,n,m

      maxA_ = -1d99
      minA_ =  1d99
      do n = 1, mat_Ads
         maxA_(n) = max(maxA_(n),maxval(mat_A(:,:,:,n)))
         minA_(n) = min(minA_(n),minval(mat_A(:,:,:,n)))
      end do
      do n = 1, mat_Ads
         maxA(n) = max(maxA(n),maxA_(n))
         minA(n) = min(minA(n),minA_(n))
      end do

      if(fls==2) then
         print*, 'saving mat_vec'//cnum//'.cdf ...'
         e=nf90_create('mat_vec'//cnum//'.cdf', nf90_clobber, f)
      else
         print*, 'saving mat_vec.cdf ... ' 
         e=nf90_create('mat_vec.cdf', nf90_clobber, f)
      end if
   
      e=nf90_def_dim(f,'xd',mat_nx,  xd)	! position dimensions
      e=nf90_def_dim(f,'yd',mat_ny,  yd)
      e=nf90_def_dim(f,'zd',mat_nz,  zd)
      e=nf90_def_dim(f,'Ad',mat_Ads, Ad)
      e=nf90_def_var(f,'x',nf90_double,(/xd/), xid)      
      e=nf90_def_var(f,'y',nf90_double,(/yd/), yid)      
      e=nf90_def_var(f,'z',nf90_double,(/zd/), zid)      
      dims = (/yd,xd,zd,Ad/)
      e=nf90_def_var(f,'A',nf90_double,dims, Aid)
      e=nf90_enddef(f)

      e=nf90_put_var(f,xid, mat_x)
      e=nf90_put_var(f,yid, mat_y)
      e=nf90_put_var(f,zid, mat_z)
      do i = 1, mat_Ads
         do n = 1, mat_nz, 100
            m = min(mat_nz,n+99)
            e=nf90_put_var(f,Aid, mat_A(:,:,n:m,i),start=(/1,1,n,i/)) 
         end do
      end do
      e=nf90_close(f)
   
   end subroutine dx_dump


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   subroutine dx_deletedata()
      deallocate(mat_x)
      deallocate(mat_y)
      deallocate(mat_z)
      deallocate(mat_A)
   end subroutine dx_deletedata
   

!*************************************************************************
 end module dx
!*************************************************************************

!*************************************************************************
#include "../parallel.h"
 PROGRAM ASH2DX
!*************************************************************************
   use io
   use dx
   implicit none
   character(200) :: path
   integer :: n, pthlen
   integer :: e, f, i
   double precision :: d
   
 9 format(A)

   print*, 'initialising...'
   call mpi_precompute()
   call par_precompute()
   call mes_precompute()
   call var_precompute()
   call tra_precompute()
   call tim_precompute()
   call vel_precompute()
   call  io_precompute()
   if(mpi_sze/=1) stop 'set _Np=1'

   call  dx_precompute()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   print*, 'Enter output variable:'
   print*, ' 1, u;'
   print*, ' 2, curl(u);'
   print*, ' 3, z . curl(u)  =  -Lap_h T = -T1'
   print*, ' 4, z . curl curl(u)  =  Lap Lap_h P = P2'
   print*, ' 5, <P2>_z  -->  h'
   print*, ' 6, P2 - h  -->  P2^'
   print*, ' 7, pressure'
   print*, ' 8, lambda_2 -- vortex'

   read(*,*) var

   print*, 'Substract mean profile?'
   print*, ' 0, use total flow;  1, sub HPf;  2, set 00mode to zero;'   
   print*, ' 3, sub from meanprof file.'
   read(*,*) submean
   if(submean==3) then
      print*, 'Enter filename (must have no header and same i_N):'
      read(*,9) path
      open(10,status='old',file=path)
      do n = 1, i_N
         read(10,*) d, velprof(n)
      end do
      close(10)
   end if

   print*, 'Shift z-frame?'
   print*, ' 0, false;  1, detect puff frame;  2, z = a t + b'
   print*, ' 3, read sz.dat file'
   read(*,*) zshft
   if(zshft==2) then
      print*, ' Enter a, b:'
      read(*,*) zshfta, zshftb
   end if
   if(zshft==1 .or. zshft==2 .or. zshft==3) then
      print*, ' Enter shift data right zoffset:'
      read(*,*) zoffset
   end if
   if(zshft==3) open(77,status='old',file='sz.dat')
   
   print*, ' Enter shift data theta offset:'
   read(*,*) toffset

   print*, 'Enter output data domain:'
   print*, '1, volume;  2, (r,th) cross-section;  3, (r,z) cross-section;'
   print*, '4, (th,z) cylindrial section.'
   read(*,*) ddom

   print*, 'i_N  =', i_N
   print*, 'i_Th =', i_Th
   print*, 'i_Z  =', i_Z
   print*, 'L    =', 2d0*d_PI/d_alpha

   if(ddom==1) then
      print*, 'Enter range z1 to z2, each in range [0,L]:'
      read(*,*) zval, zval2
      print*, 'Enter nx, ny, nz:'
      read(*,*) mat_nx, mat_ny, mat_nz      
   else if(ddom==2) then
      print*, 'Enter z value in range [0,L]:'
      read(*,*) zval
      print*, 'Enter nx, ny:'
      read(*,*) mat_nx, mat_ny
      zval2  = zval
      mat_nz = 1
   else if(ddom==3) then
      print*, 'Enter th value in range [0,1] (mapped to [0,2pi]):'
      read(*,*) thval
      thval = thval * 2d0*d_PI
   else if(ddom==4) then
      print*, 'Enter r value in range [0,1]:'
      read(*,*) rval
   end if

   if(var==1 .or. var==2) then
      mat_Ads = 3
      print*, 'Enter ouput vectors:'
      print*, '1, (Ax,Ay,Az);  2, (Ar,Ath,Az).'
      read(*,*) dvec
   else
      mat_Ads = 1
      dvec    = 2
   end if

   print*, 'Multiple files?:  1, single file;  2, series'
   read(*,*) fls
   if(fls==1) then
      print*, 'Enter state file:'
      read(*,9) io_statefile   
      ffrst = 1 
      flast = 1
      fstep = 1
   else 
      fls = 2
      print*, 'Enter path to files:'
      read(*,9) path
      pthlen = 1
      do while(path(pthlen+1:pthlen+1)/=' ')
         pthlen = pthlen+1
      end do
      print*, 'Enter file numbers first,last,step:'
      read(*,*) ffrst, flast, fstep
   end if

   tim_dt = 0.01d0
   call vel_matrices()
   
   maxA = -1d99
   minA =  1d99

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
   open(11,status='unknown',file='mat_maxmin')

   do fl = ffrst, flast, fstep
      
      if(fls==2) then
         write(cnum,'(I4.4)') fl
         io_statefile = path(1:pthlen)//'state'//cnum//'.cdf.dat'   
         print*,'loading '//path(1:pthlen)//'state'//cnum//'.cdf.dat'
      end if
      dtprev = tim_dt
      call io_load_state()
      print*, 'processing...'
      print*, ' fl = ', fl
      
      call dx_preparedata()
         
      if(ddom==1 .or. ddom==2) then
         call dx_volume()
      else if(ddom==3) then
         call dx_rz_crsect()
      else if(ddom==4) then
         call dx_thz_cylsect()
      end if

      call dx_dump()
      call dx_deletedata()

      if(fls==2)  write(11,'(I4,6e16.8)') fl,  &
         (maxA_(n), n=1,mat_Ads), (minA_(n), n=1,mat_Ads)

   end do

   write(*,*) 'writing to mat_maxmin:'
   write(*,*) 'max: ', (maxA(n), n=1,mat_Ads)
   write(*,*) 'min: ', (minA(n), n=1,mat_Ads)
   write(11,*) '# data max:', (maxA(n), n=1,mat_Ads)
   write(11,*) '# data min:', (minA(n), n=1,mat_Ads)   
   write(11,*) '# d_alpha : ',d_alpha
   write(11,*) '# i_Mp : ',i_Mp
   write(11,*) '# i_N : ', i_N
   write(11,*) '# i_Th: ', i_Th
   write(11,*) '# i_Z : ', i_Z
   close(11)

!*************************************************************************
 END PROGRAM ASH2DX
!*************************************************************************

!----------------------------------------
! on exit E_ eigvals in descending order
!----------------------------------------
 subroutine get_eigvals(N,A,LDA,E_)
    integer, intent(in) :: N, LDA
    double precision, intent(inout) :: A(LDA,N), E_(N)
    double precision :: d(N), e(N-1), tau(N-1), work(N*N)
    integer :: info
    
    call dsytrd('U',N,A,LDA,d,e,tau,work,N*N,info)
    if(info/=0) stop 'dsytrd'
    call dsterf(N,d,e,info)
    if(info/=0) stop 'dsterf'
    
    E_ = d(N:1:-1)    
 
 end subroutine get_eigvals



