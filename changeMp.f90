!*************************************************************************
 PROGRAM CHANGEMP
!*************************************************************************
   use netcdf
   implicit none
   character(200)   :: fname
   integer :: e,f,i
   integer :: Mp,m
 9 format(A)
   
   print*, 'Enter state filename (WARNING: OVERWRITES FILE):'
   read(*,9) fname

   e=nf90_open(fname, nf90_write, f)
   if(e/=nf90_noerr) stop 'file not found!'
   
   print*, 'Output Mp'
   read(*,*) Mp

   e=nf90_inq_varid(f,'Ur', i)
   e=nf90_get_att(f,i,'Mp', m)
   e=nf90_put_att(f,i,'Mp', Mp)
   e=nf90_inq_varid(f,'Ut', i)
   e=nf90_put_att(f,i,'Mp', Mp)
   e=nf90_inq_varid(f,'Uz', i)
   e=nf90_put_att(f,i,'Mp', Mp)

   print*, 'Mp: ', m, ' --> ', Mp
   
   e=nf90_close(f)

!*************************************************************************
 END PROGRAM CHANGEMP
!*************************************************************************
