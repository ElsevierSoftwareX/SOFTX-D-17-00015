function [] = isosurf(ncfile,afac)

% crosssec can be called without any arguments:
if nargin < 2 ; afac = 0.8 ; end
if nargin < 1 ; ncfile = 'mat_vec.cdf'; end
close all;
figure('Position',[0 0 880 880]);


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
size(z)
% A has size (ny,nx,nz,Ads)
% Ads=1 for scalar data, =3 for vector data


[X,Y,Z] = meshgrid(x,y,z);

amin = min(min(min(A(:,:,:,1))))
amax = max(max(max(A(:,:,:,1))))
aiso = (1.-afac)*amin + amax*afac

p = patch(isosurface(X,Y,Z,A(:,:,:,1),aiso));
isonormals(X,Y,Z,A(:,:,:,1),p)
set(p,'FaceColor','blue','EdgeColor','none');
daspect([1,1,1])
axis tight;
axis([-1 1 -1 1]);
hold on;

zlim = get(gca,'ZLim');
theta=0:0.01:2*pi;
points(1,:) = cos(theta);
points(2,:) = sin(theta);
points(3,:) = theta .* 0 + zlim(1);
points(4,:) = theta .* 0 + zlim(2);
plot3(points(1,:),points(2,:),points(3,:),'r-');
plot3(points(1,:),points(2,:),points(4,:),'r-');

camlight 
lighting gouraud
%view(3);
%view([0 90]); camproj('orthographic'); camlight(-25,45);
view([160 -50]); camproj('perspective'); camroll(270);


end
