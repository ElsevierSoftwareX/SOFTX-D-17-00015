#ifndef parallel_h
#define parallel_h

/***************************************************************/
#define _Nr 1
#define _Ns 1
#define _CPUTIME
/***************************************************************/

/*-------------------------------------------------------------*/
#define _Np (_Nr*_Ns)
#if _Np != (1*1)
#define _MPI
#endif

#define _Ms  (i_M/_Ns)
#define _Ms1 (_Ms-1)
#define _Hs1 ((2*i_K-1)*_Ms-1)

/*-------------------------------------------------------------*/
#if _Np != 1
#define _MPI
#endif

/*-------------------------------------------------------------*/
#ifdef _MPI
/*-------------------------------------------------------------*/

#define _loop_km_vars  integer :: pH0,pH1,nh,k,m

#define _loop_km_begin  /*
*/ pH0=var_H%pH0;  pH1=var_H%pH1;  /*
*/ m = (pH0+i_K1)/(2*i_K-1);  /*
*/ k = pH0 - m*(2*i_K-1) - 1;  /*
*/ do nh = 0, pH1;  /*
*/  k = k + 1;  /*
*/  if(k==i_K) then;  /*
*/   k = -i_K1;  m = m + 1;  /*
*/  end if

#define _loop_km_end   end do

/*-------------------------------------------------------------*/
#else  /* ndef _MPI */
/*-------------------------------------------------------------*/

#define _loop_km_vars  integer :: nh,k,m

#define _loop_km_begin /*
*/ nh = -1;  /*
*/ do m = 0, i_M1;  /*
*/  do k = -i_K1, i_K1;  /*
*/   if(k<0 .and. m==0) cycle;  /*
*/   nh = nh  + 1
         
#define _loop_km_end  end do;  end do

/*-------------------------------------------------------------*/
#endif /* def _MPI */
/*-------------------------------------------------------------*/

#endif /* ndef parallel_h */
