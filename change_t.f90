!*************************************************************************
 PROGRAM CHANGET
!*************************************************************************
   use netcdf
   implicit none
   character(200)   :: fname
   double precision :: d
   integer :: e,f,i
 9 format(A)
   
   print*, 'Enter state filename: (WARNING -- OVERWRITES FILE)'
   read(*,9) fname

   e=nf90_open(fname, nf90_write, f)
   if(e/=nf90_noerr) stop 'file not found!'
   
   e=nf90_get_att(f,nf90_global,'t', d)
   
   print*, 'Enter new t:     ( was,',real(d),')'
   read(*,*) d   
   e=nf90_put_att(f,nf90_global,'t', d)
   
   e=nf90_close(f)
   
!*************************************************************************
 END PROGRAM CHANGET
!*************************************************************************
