function [r] = getGRID(Resolution,Lx,Ly,Lz,Rc)
%%

% Resolution = Diameter/nX;
dx = Resolution;
dy = dx;
dz = dx;

Nx = ceil(Lx/dx);
Ny = ceil(Ly/dy);
Nz = ceil(Lz/dz);

%
xmax = (Nx/2)*dx - dx/2;
ymax = (Ny/2)*dy - dy/2;
zmax = (Nz/2)*dz - dz/2;

%
x = -xmax:dx:xmax;
y = -ymax:dy:ymax;
z = -zmax:dz:zmax;

% Move the center from origin
if nargin > 4
    x = x + Rc(1);
    y = y + Rc(2);
    z = z + Rc(3);
end

%
r = grid3d(x,y,z);

end

function r = grid3d(x, y, z)
%%

if nargin == 1
    y = x;
    z = x;
end


L = length(x);
M = length(y);
N = length(z);

r = zeros(L,M,N,3);

for ix = 1:L
    xx = x(ix);
    for iy = 1:M
        yy = y(iy);
        for iz = 1:N
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end

end
