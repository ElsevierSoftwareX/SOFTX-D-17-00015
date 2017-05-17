!*************************************************************************
!  Compares (Mp,alpha) of parameters.f90 and of input file to determine
!  either number of copies to make or extraction of modes.
!  Note that N,K,M of parameters.f90 must be sufficient to avoid 
!  trunctaion of either input or output: compare K-1 with (K-1)*Kp etc.
!*************************************************************************
#include "../parallel.h"
 PROGRAM MAIN
!*************************************************************************
   use io
   implicit none
   integer :: e,f,i, ex
   integer :: n, Kp_,Mp_, Mp__
   type (coll) :: c1,c2,c3
   double precision :: alp
   _loop_km_vars
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
   if(mpi_sze/=1)  stop 'Np_'
   
   print*, 'Output parameters:'
   print*, 'K,M   = ', i_K,i_M
   print*, 'Mp    = ', i_Mp
   print*, 'alpha = ', d_alpha
   
   print*, 'Enter state filename'
   print*, '(WARNING: writes to state0000.cdf.dat, vel_spec0000.dat):'
   read(*,9) io_statefile
   e=nf90_open(io_statefile, nf90_nowrite, f)
   if(e/=nf90_noerr) stop 'file not found!'
   e=nf90_get_att(f,nf90_global,'alpha', alp)
   e=nf90_inq_varid(f,'Ur', i)
   e=nf90_get_att(f,i, 'K', k)
   e=nf90_get_att(f,i, 'M', m)
   e=nf90_get_att(f,i,'Mp', Mp__)
   e=nf90_close(f)

   print*, 'Input parameters:'
   print*, 'K,M   = ', k,m
   print*, 'Mp    = ', Mp__
   print*, 'alpha = ', alp
   
   Kp_ = int(max(alp/d_alpha,d_alpha/alp)+0.5d0)
   Mp_ = int(max(dble(Mp__)/dble(i_Mp),dble(i_Mp)/dble(Mp__))+0.5d0)
   ex = 1
   if(d_alpha>alp*1.1d0 .or. i_Mp>Mp__) ex = 0
   print*, 'Mp: ', Mp_, ' copies/extract theta' 
   print*, 'Kp: ', Kp_, ' copies/extract z' 

   call system('cp '//io_statefile//' TEMP.cdf')
   e=nf90_open('TEMP.cdf', nf90_write, f)
   if(e/=nf90_noerr) stop 'TEMP file not found!'
   e=nf90_inq_varid(f,'Ur', i)
   e=nf90_put_att(f,i,'Mp', i_Mp)
   e=nf90_inq_varid(f,'Ut', i)
   e=nf90_put_att(f,i,'Mp', i_Mp)
   e=nf90_inq_varid(f,'Uz', i)
   e=nf90_put_att(f,i,'Mp', i_Mp)
   e=nf90_close(f)
   io_statefile = 'TEMP.cdf'
   if(mpi_rnk==0)  print*, 'loading state...'
   call io_load_state()
   call system('rm TEMP.cdf')

   tim_t = 0d0
   call var_coll_copy(vel_ur, c1)
   call var_coll_copy(vel_ut, c2)
   call var_coll_copy(vel_uz, c3)
   call var_coll_init(vel_ur)
   call var_coll_init(vel_ut)
   call var_coll_init(vel_uz)
   
   _loop_km_begin
      if(    m*Mp_ >i_M1) cycle
      if(abs(k*Kp_)>i_K1) cycle
      call harmonic(k*Kp_,m*Mp_, n)
      if(ex==1) then
         vel_ur%Re(:,n) = c1%Re(:,nh)
         vel_ur%Im(:,n) = c1%Im(:,nh)
         vel_ut%Re(:,n) = c2%Re(:,nh)
         vel_ut%Im(:,n) = c2%Im(:,nh)
         vel_uz%Re(:,n) = c3%Re(:,nh)
         vel_uz%Im(:,n) = c3%Im(:,nh)
      else
         vel_ur%Re(:,nh) = c1%Re(:,n)
         vel_ur%Im(:,nh) = c1%Im(:,n)
         vel_ut%Re(:,nh) = c2%Re(:,n)
         vel_ut%Im(:,nh) = c2%Im(:,n)
         vel_uz%Re(:,nh) = c3%Re(:,n)
         vel_uz%Im(:,nh) = c3%Im(:,n)
      end if
   _loop_km_end
   
   print*, ' saving vel_spec0000.dat'
   call io_save_spectrum()
   call io_save_state()

         
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

 subroutine harmonic(k,m, nh)
    use parameters
    implicit none
    integer, intent(in)  :: k,m
    integer, intent(out) :: nh

    nh = m*(2*i_K1+1) + k

 end subroutine harmonic

