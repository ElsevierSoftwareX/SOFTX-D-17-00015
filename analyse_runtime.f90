!*************************************************************************
 module runtimedata
   use variables
   implicit none
   save
   integer :: p1, p2
   double precision :: d1, d2, d3
   type (phys) :: F
 end module runtimedata


!*************************************************************************
#define var_null util
#include "../program/main.f90"
!*************************************************************************
 subroutine util(FLAG)
   use runtimedata
   implicit none
   integer, intent(in) :: FLAG
   logical, save :: set = .false.
				! this subroutine is called every timestep, 
				! just before writing usual data.
   if(.not.set) then
      set = .true.
				! initialise
      open(99,status='old',file='myparams.in')
      read(99,*) p1, p2
      close(99)      
   
      open(61,status='unknown',file='mydata.dat')
   end if

   if(FLAG==1) then
                                ! called just before writing data in main
      ! analyse state,
      ! write data to unit 61...

   else if(FLAG==2) then
                                ! called after eval of nonlin terms
      ! add extra terms to equation,
      ! use user parameters p1, p2...
      
   end if


 end subroutine util
