function [x,r,th,x_1D,r_1D,nx,nr,nth] = read_CylGrid_656_138_128_xyz(file)

nx      = 656;
ny      = 138;
nz      = 128;

nr      = ny;
nth     = nz;

fileID  = fopen(file);
np      = fread(fileID,1,'int');
xyz     = fread(fileID,[3,np],'float');
fclose(fileID);

x       = reshape(squeeze(xyz(1,:)),nz,ny,nx);
y       = reshape(squeeze(xyz(2,:)),nz,ny,nx);
z       = reshape(squeeze(xyz(3,:)),nz,ny,nx);

x_1D    = squeeze(x(1,1,:));

th      = 2*pi*[0:nth-1]'/nth;
r_1D    = squeeze(y(1,:,1));
[x r]   = meshgrid(x_1D,r_1D);

end
