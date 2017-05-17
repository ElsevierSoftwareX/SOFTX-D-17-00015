function [] = crosssec(nlev,ncfile)

% crosssec can be called without any arguments:
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

[X,Y] = meshgrid(x,y);

amin = min(min(min(A(:,:,1,3))))
amax = max(max(max(A(:,:,1,3))))
step = (amax-amin)/nlev ;
levels = amin:step:amax ;

contourf(X,Y,A(:,:,1,3),levels);
%contour(X,Y,A(:,:,1,3),levels);
colormap('hot');
brighten(0.4);

hold on;
xx = x(1:2:end);
yy = y(1:2:end);
Ax = A(1:2:end,1:2:end,1,1);
Ay = A(1:2:end,1:2:end,1,2);
quiver(xx,yy,2*Ax,2*Ay,0,'color','b');

if x(end) < 1.
   theta=0:0.001:2*pi;
   xxx=cos(theta);
   yyy=sin(theta);
   plot(xxx,yyy);  
   axis square;
end

hold off;


end
