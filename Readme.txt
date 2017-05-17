-----------------------------------------------------------------------
The utility  utils/prim2matlab.f90  
creates NetCDF (.cdf) files, which are readable by matlab, visit, etc.

Included in this directory, matlab/ 
are a few matlab scripts that use output from prim2matlab.out

-----------------------------------------------------------------------
Usage:

In parallel.h set _Nr and _Ns both to 1 .
In Makefile change UTIL to prim2matlab,
ensure using serial COMPILER.
   make
   make util

Move prim2matlab.out to your run directory.

Run prim2matlab.out, enter your options:
   The data may be shifted/rotated before output;
      enter 0. for no offset.
   The data must be interpolated onto a Cartesian
      grid; pick number of pts e.g. for nx, ny
      50 50 
   Given the data is on a Cartesian grid, usally
      want Cartesian vector data.

After choosing options, prim2matlab creates 
   mat_maxmin,
an info file, useful for, e.g. selecting isosurface 
values, and a data file for matlab,
   mat_vec.cdf .

With mat_vec.cdf the in current directory, e.g.,
   cp ~/svn/pipe_primitive/matlab/crosssec.m .
   matlab &
   > crosssec
or
   > crosssec(20)		for 20 isocontours
   > crosssec(20,'file.cdf') 	for another file

----------------------------------------------------------------------
Format of mat_vec.cdf:

NetCDF format with variables A (data), x,y,z (grid points).
   x has dimension x(nx)
   A has dimension A(nx,ny,nz,Ads)

nz=1  for a cross-section.
Ads=1 for scalar data.
Ads=3 for a vector.

