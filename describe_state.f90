!*************************************************************************
 PROGRAM DESC
!*************************************************************************
   use netcdf
   implicit none
   character(200)   :: fname
   integer :: N,K,M,Mp
   integer :: e,f,i,i_
   double precision :: t,alpha,Re
 9 format(A)
   
   print*, 'Enter state filename:'
   read(*,9) fname

   e=nf90_open(fname, nf90_nowrite, f)
   if(e/=nf90_noerr) stop 'file not found!'
   
   e=nf90_get_att(f,nf90_global,'t', t)
   print*, 't     = ', t
   e=nf90_get_att(f,nf90_global,'Re', Re)
   print*, 'Re    = ', Re
   e=nf90_get_att(f,nf90_global,'alpha', alpha)
   print*, 'alpha = ', alpha

   e=nf90_inq_dimid(f, 'r', i)
   e=nf90_inquire_dimension(f, i, len=N)

   e=nf90_inq_varid(f,'P2', i_)
   if(e==nf90_noerr) then
      i = i_
      print*, 'Tor-Pol format'
      e=nf90_inq_varid(f,'T', i_)
      print*, '   save_all = ', (e==nf90_noerr)
   end if
   e=nf90_inq_varid(f,'Ur', i_)
   if(e==nf90_noerr) then
      i = i_
      print*, 'Primitive variable format'
   end if
   e=nf90_get_att(f,i,'K',  K)
   e=nf90_get_att(f,i,'M',  M)
   e=nf90_get_att(f,i,'Mp', Mp)

   print*, 'i_N   = ', N
   print*, 'i_K   = ', K
   print*, 'i_M   = ', M
   print*, 'i_Mp  = ', Mp
   
   e=nf90_close(f)
   

!*************************************************************************
 END PROGRAM DESC
!*************************************************************************
