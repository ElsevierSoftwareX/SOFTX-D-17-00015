function [] = streamfn(nlev,ncfile)

if nargin < 2 ; ncfile = 'mat_vec.cdf'; end
if nargin < 1 ; nlev = 10 ; end
close all;


fid = netcdf.open(ncfile,'NC_NOWRITE');
xid = netcdf.inqVarID(fid,'x');
yid = netcdf.inqVarID(fid,'y');
zid = netcdf.inqVarID(fid,'z');
Aid = netcdf.inqVarID(fid,'A');

x = netcdf.getVar(fid,xid);
y = netcdf.getVar(fid,yid);
z = netcdf.getVar(fid,zid);
A = netcdf.getVar(fid,Aid);
netcdf.close(fid);
size(x)
size(y)
z
% z is the value at which the cross section was taken.
% A has size (ny,nx,1,Ads)
% Ads=1 for scalar data, =3 for vector data.
% Expecting Ads=1, scalar streamfunction.

[X,Y] = meshgrid(x,y);

% amin = 0. 
amin = min(min(min(A(:,:,1,1))))
amax = max(max(max(A(:,:,1,1))))
%step = (amax-amin)/nlev ;
%levels = amin:step:amax ;
nstep = (0.001-amin)/(nlev);
pstep = (amax-0.001)/(nlev/4.);
nlevels = amin:nstep:0.001 ;
plevels = 0.001:pstep:amax ;
levels = cat(2, nlevels, plevels) ;

contourf(X,Y,A(:,:,1,1),levels);
% contour(X,Y,A(:,:,1,1),levels);
colormap('hot');
brighten(0.4);

hold on;

hold off;

end
