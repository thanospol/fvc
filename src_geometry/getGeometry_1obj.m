function [r,EMT] = getGeometry_1obj(points_min,ObjProperties)
% Discretization for various geometry
%    
% -------------------------------------------------------------------------
% Basic geometries with homogeneous temperature and permittivity distribution:
% shape = 'Sphere': a sphere with R=Radius.
%         'Cube':  a cube with edge-length 2*Radius.
%         'Ellipsoid': a ellipsoid with semi-axis (major in x,y,z directions): Radius*SIZE(1,2,3)
%         'Cylinder': a cylinder, with ellipse x/y radius and half-height to be Radius*SIZE(1,2,3) respectivley.
%    
% -------------------------------------------------------------------------
% Basic geometries with linearly varying temperature and permittivity along z-direction:
%         T=Temp(1)+(Temp(2)-Temp(1))*(z+R)/2R;Er=e_r(1)+(e_r(2)-e_r(1))*(z+R)/2R    
% shape = 'Cubelinearz': for a cube
%         'Spherelinearz': for a sphere
%         'Ellispoidlinearz': for a ellipsoid 
%    
% -------------------------------------------------------------------------    
% Basic geometries with linearly varying temperature and permittivity along radius-direction:
%         T=Temp(1)+(Temp(2)-Temp(1))*r/R;Er=e_r(1)+(e_r(2)-e_r(1))*r/R        
% shape = 'SpherelinearR': for a sphere
%    
% -------------------------------------------------------------------------
% Complex geometries with complex temperature/permittivity distributions:
% shape = 'Hemispheroidshell': a hemispheroid with shell. The core: Temp(1)/e_r(1); the shell: Temp(2)/e_r(2)
%                              The size for the core: semi-axis (in x,y,z directions) Radius*SIZE(1,2,3)
%                              the shell: a relative thickness SIZE(4)
%         'Twosphere': two spheres with the same r=Radius, and the ditance between the nearest point Radius*SIZE(1).
%                      The left/right sphere: Temp(1/2), and e_r(1/2).
%         'Mushroom':  a mushroom-shape:the spheroid: Temp(1)/e_r(1); the cylinder: Temp(2)/e_r(2)
%                     The size for the spheroid: semi-axis length (in x,y,z directions) Radius*SIZE(1,2,3)
%                     The cylinder: radius Radius*SIZE(4), and height Radius*SIZE(5)
%         'Ellipsoid_mix2':  an ellipsoid has Temp(1)/e_r(1) for the part Z>Radiusz-2*Radiusz*SIZE(4), and Temp(2)/e_r(2) for 
%                            the rest.  Three semi-axis lengths are (major in x,y,z directions): Radius*SIZE(1,2,3)
%         'Hemispheroidshell_Tfromdata': a hemispheroid with shell: the temperature profile is from comsol simulation data.
%                             The permittivity of the core is a function of temperature e_r=core_diperseT(T).
%                             The size: core has semi-axis SIZE(1,2,3) for x,y,z directions,
%                             and the shell has the relative thickness SIZE(4)
%                             The permittivity of the  shell is e_r
    
Epsilon      = ObjProperties.Epsilon;
Mu           = ObjProperties.Mu;
shape        = ObjProperties.shape;
Radius       = ObjProperties.Radius;
if isfield(ObjProperties,'SIZE')
    SIZE=ObjProperties.SIZE;
end
if isfield(ObjProperties,'Temp_profile')
    Temp_profile = ObjProperties.Temp_profile;
end
% -------------------------------------------------------------------------
% DEFINE GEOMETRY
% -------------------------------------------------------------------------
if strcmp(shape,'Cylinder')||strcmp(shape,'Ellipsoid')||strcmp(shape,'Ellipsoidlinearz')||strcmp(shape,'Ellipsoid_mix2')

    Rc = [0; 0; 0];
    
    Diam_x = 2*Radius*SIZE(1);  % Ellipse x diameter
    Diam_y = 2*Radius*SIZE(2);  % Ellipse y diameter
    Diam_z = 2*Radius*SIZE(3);  % Height
    
    Diameter_min = min([Diam_x,Diam_y,Diam_z]);
    
    Lx = Diam_x;
    Ly = Diam_y;
    Lz = Diam_z;
end

if strcmp(shape,'Sphere')||strcmp(shape,'Cube')||strcmp(shape,'Cubelinearz')||strcmp(shape,'Spherelinearz')...
        ||strcmp(shape,'SpherelinearR')
    Rc = [0; 0; 0];
    
    Diameter_min = 2*Radius;
    
    Lx = Diameter_min;
    Ly = Diameter_min;
    Lz = Diameter_min;

end

if strcmp(shape,'Hemispheroidshell')||strcmp(shape,'Hemispheroidshell_Tfromdata')

    Rc = [0; 0; 0];

    Diam_x = 2*Radius*SIZE(1)*(1+SIZE(4));
    Diam_y = 2*Radius*SIZE(2)*(1+SIZE(4));
    Diam_z = 2*Radius*SIZE(3)*(1+SIZE(4));

    Diameter_min = min([Diam_x,Diam_y,Diam_z]);

    Lx = Diam_x;
    Ly = Diam_y;
    Lz = Diam_z;
end

if strcmp(shape,'Twosphere')

    Rc = [0; 0; 0];
    
    Diam_x = Radius*(4*SIZE(1)+SIZE(2));
    Diam_y = 2*Radius*SIZE(1);
    Diam_z = 2*Radius*SIZE(1);
    
    Diameter_min = min([Diam_x,Diam_y,Diam_z]);
    
    Lx = Diam_x;
    Ly = Diam_y;
    Lz = Diam_z;
end        

if strcmp(shape,'Mushroom')
    
    Rc = [0; 0; 0];

    Diam_x = 2*Radius*max(SIZE(1),SIZE(4));
    Diam_y = 2*Radius*SIZE(2);
    Diam_z = 2*Radius*max(SIZE(3),SIZE(5));

    Diameter_min = min([Diam_x,Diam_y,Diam_z]);

    Lx = Diam_x;
    Ly = Diam_y;
    Lz = Diam_z;
end
% -------------------------------------------------------------------------
% DISCRETIZATION
% -------------------------------------------------------------------------
Resolution = Diameter_min/points_min;

[r] = getGRID(Resolution,Lx,Ly,Lz,Rc);

if strcmp(shape,'Sphere')
    [Er,Mr,Tr,idx] = sphere(r,Rc,Radius,Epsilon,Mu,Temp_profile);
end
if strcmp(shape,'Spherelinearz')
    [Er,Mr,Tr,idx] = sphere(r,Rc,Radius,Epsilon,Mu,Temp_profile);
    Z=r(:,:,:,3)-Rc(3,1);
    Er(idx) = Epsilon(1)+ (Z(idx)+Radius)*(Epsilon(2)-Epsilon(1))/(2*Radius);
    Tr(idx) = Temp_profile(1)+ (Z(idx)+Radius)*(Temp_profile(2)-Temp_profile(1))/(2*Radius);
end
if strcmp(shape,'SpherelinearR')
    [Er,Mr,Tr,idx] = sphere(r,Rc,Radius,Epsilon,Mu,Temp_profile);
    radius_dis=sqrt((r(:,:,:,1)-Rc(1,1)).^2+(r(:,:,:,2)-Rc(2,1)).^2+(r(:,:,:,3)-Rc(3,1)).^2);
    Er(idx) = Epsilon(1) + (Epsilon(2)-Epsilon(1))*radius_dis(idx)/Radius;
    Tr(idx) = Temp_profile(1) + (Temp_profile(2)-Temp_profile(1))*radius_dis(idx)/Radius;
end


if strcmp(shape,'Cube')
    [Er,Mr,Tr] = cube(r,Rc,Radius,Epsilon,Mu,Temp_profile);
end
if strcmp(shape,'Cubelinearz')
    [Er,Mr,Tr,idx] = cube(r,Rc,Radius,Epsilon,Mu,Temp_profile);
    Z=r(:,:,:,3)-Rc(3,1);
    Er(idx) = Epsilon(1)+ (Z(idx)+Radius)*(Epsilon(2)-Epsilon(1))/(2*Radius);
    Tr(idx) = Temp_profile(1)+ (Z(idx)+Radius)*(Temp_profile(2)-Temp_profile(1))/(2*Radius);
end


if strcmp(shape,'Ellipsoid')
    [Er,Mr,Tr] = ellipsoid(r,Rc, Radius,Epsilon,Mu,Temp_profile,SIZE);
end
if strcmp(shape,'Ellipsoidlinearz')
    [Er,Mr,Tr,idx] = ellipsoid(r,Rc,Radius,Epsilon,Mu,Temp_profile,SIZE);
    Z=r(:,:,:,3)-Rc(3,1);
    Er(idx) = Epsilon(1)+ (Z(idx)+Radius)*(Epsilon(2)-Epsilon(1))/(2*Radius);
    Tr(idx) = Temp_profile(1)+ (Z(idx)+Radius)*(Temp_profile(2)-Temp_profile(1))/(2*Radius);
end


if strcmp(shape,'Cylinder')
    [Er,Mr,Tr] = cylinder(r,Rc, Radius,Epsilon,Mu,Temp_profile,SIZE);
end

if strcmp(shape,'Hemispheroidshell')
    [Er,Mr,Tr] = hemispheroidshell(r,Rc,Radius,Epsilon,Mu,Temp_profile,SIZE);
end
if strcmp(shape,'Hemispheroidshell_Tfromdata')
    [Er,Mr,Tr] = hemispheroidshell_Tfromdata(r,Rc,Radius,Epsilon,Mu,SIZE,ObjProperties.dataT,ObjProperties.core_diperseT);
end

if strcmp(shape,'Twosphere')
    [Er,Mr,Tr] = twosphere(r,Rc,Radius,Epsilon,Mu,Temp_profile,SIZE);
end

if strcmp(shape,'Mushroom')
    [Er,Mr,Tr] = mushroom(r,Rc,Radius,Epsilon,Mu,Temp_profile,SIZE);
end

if strcmp(shape,'Ellipsoid_mix2')
    [Er,Mr,Tr] = ellipsoid_mix2(r,Rc,Radius,Epsilon,Mu,Temp_profile,SIZE);
end

EMT.Er = Er;
EMT.Mr = Mr;
EMT.Tr = Tr;
end
%% local functions: shapes
function [epsilon_r,mu_r,T,idx] = sphere(r, Rc, Radius, e_r, m_r, Temp)
%   Function that generates a homogeneous sphere with R=Radius

[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = zeros(L,M,N);
% define sphere
object = @(r)( (r(:,:,:,1) - Rc(1,1) ).^2 + ( r(:,:,:,2) - Rc(2,1) ).^2 + ( r(:,:,:,3) - Rc(3,1) ).^2 < Radius^2) ;
pointobject = object(r); % ones for domain elements in the sphere
idx = find(pointobject(:)); % get indexes of elements

epsilon_r(idx) = e_r(1); 
mu_r(idx) = m_r; 
T(idx) = Temp(1);
end
function [epsilon_r,mu_r,T,idx] = cube(r, Rc, Radius, e_r, m_r, Temp)
%   Function that generates a homogeneous cube with edge-length 2*Radius
[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = zeros(L,M,N);

% define cube
object = @(r)( (abs(r(:,:,:,1)-Rc(1,1))<=Radius) & (abs(r(:,:,:,2)-Rc(2,1))<=Radius) & (abs(r(:,:,:,3)-Rc(3,1))<=Radius));
pointobject = object(r); % ones for domain elements in the cube
idx = find(pointobject(:)); % get indexes of elements

epsilon_r(idx) = e_r(1); 
mu_r(idx) = m_r; 
T(idx) = Temp(1);
end
function [epsilon_r,mu_r,T,idx] = ellipsoid(r,Rc,Radius,e_r,m_r,Temp,SIZE)
% Function that generates a homogeneous ellipsoid with semi-axis (major in x,y,z directions): Radius*SIZE(1,2,3)

Radiusx=Radius*SIZE(1);
Radiusy=Radius*SIZE(2);
Radiusz=Radius*SIZE(3);

[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = zeros(L,M,N);

% define ellipsoid
object = @(r)( (r(:,:,:,1)-Rc(1,1)).^2/Radiusx^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Radiusy^2 + (r(:,:,:,3)-Rc(3,1)).^2/Radiusz^2 < 1) ;
pointobject = object(r); % ones for domain elements in the shell
idx = find(pointobject(:)); % get indexes of elements
mu_r(idx) = m_r; 
epsilon_r(idx) = e_r(1);
T(idx) = Temp(1);
end
function [epsilon_r,mu_r,T,idx] = cylinder(r, Rc, Radius, e_r, m_r,Temp,SIZE)
% Function that generates a homogeneous cylinder, with ellipse radius (x,y) and half height to be Radius*SIZE(1,2,3)
    
a=Radius*SIZE(1);
b=Radius*SIZE(2);
h_half=Radius*SIZE(3);

[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = zeros(L,M,N);

% define cylinder
object = @(r)( (abs( (r(:,:,:,1)-Rc(1,1))/a + 1i*(r(:,:,:,2)-Rc(2,1))/b) < 1) & (abs((r(:,:,:,3)-Rc(3,1))) < h_half) );

pointobject = object(r); % ones for domain elements in the cylinder
idx = find(pointobject(:)); % get indexes of elements

epsilon_r(idx) = e_r(1); 
mu_r(idx) = m_r; 
T(idx) = Temp(1);
end
function [epsilon_r,mu_r,T,idx] = twosphere(r,Rc,Radius,e_r,m_r,Temp,SIZE)
%   Function that generates two spheres, with the same radius
%  , and the ditance between the nearest point SIZE(1).
%   The left/right sphere has T/E(1,2)

distance=Radius*SIZE(1);
r0=Radius;
x1=-distance/2-r0; 
x2=distance/2+r0; % center(x) for two spheres.

[L,M,N,~] = size(r);
epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N)*m_r;
T = zeros(L,M,N);

% define sphere1
sphere1 = @(r)( (r(:,:,:,1)-x1-Rc(1,1)).^2/r0^2 +  (r(:,:,:,2)-Rc(2,1)).^2/r0^2 + (r(:,:,:,3)-Rc(3,1)).^2/r0^2 < 1) ;
pointsphere = sphere1(r); 
idx = find(pointsphere(:)); % get indexes of elements
epsilon_r(idx) = e_r(1); 
T(idx) = Temp(1);

% define sphere2
sphere2 = @(r)( (r(:,:,:,1)-x2-Rc(1,1)).^2/r0^2 +  (r(:,:,:,2)-Rc(2,1)).^2/r0^2 + (r(:,:,:,3)-Rc(3,1)).^2/r0^2 < 1) ;
pointsphere = sphere2(r);
idx = find(pointsphere(:)); % get indexes of elements
epsilon_r(idx) = e_r(2); 
T(idx) = Temp(2);
end
function [epsilon_r,mu_r,T,idx] = hemispheroidshell(r,Rc,Radius,e_r,m_r,Temp,SIZE)
%   Function that generates a hemispheroid with shell.

%   The size: core has semi-axis SIZE(1,2,3) for x,y,z directions,
%   and the shell has the relative thickness SIZE(4)

%   The core has the temperature/permittivity T/e_r(1); the shell T/e_r(2)

[L,M,N,~] = size(r);

Rx=Radius*SIZE(1);
Ry=Radius*SIZE(2);
Rz=Radius*SIZE(3); %semi-axis for the core

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N)*m_r;
T = zeros(L,M,N);

% define core
core = @(r)( ((r(:,:,:,1)-Rc(1,1)).^2/Rx^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Ry^2 + (r(:,:,:,3)-Rc(3,1)).^2/Rz^2 < 1) & ((r(:,:,:,3)-Rc(3,1))>0));
pointcore = core(r);
idx = find(pointcore(:)); % get indexes of elements
epsilon_r(idx)=e_r(1);
T(idx)=Temp(1);

% define shell
shell = @(r)( ((r(:,:,:,1)-Rc(1,1)).^2/Rx^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Ry^2 + (r(:,:,:,3)-Rc(3,1)).^2/Rz^2 <= (1+SIZE(4))^2)... 
              &((r(:,:,:,1)-Rc(1,1)).^2/Rx^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Ry^2 + (r(:,:,:,3)-Rc(3,1)).^2/Rz^2 >= 1) & ...
              ((r(:,:,:,3)-Rc(3,1))>0) );
point=shell(r);
idx=find(point(:));
epsilon_r(idx) = e_r(2); 
T(idx) = Temp(2);
end
function [epsilon_r,mu_r,T,idx] = mushroom(r,Rc,Radius,e_r,m_r,Temp,SIZE)
% Function that generates a mushroom-shape, with the head a spheroid, and tail a cylinder.
% The spheroid has its semi-axis length (in x,y,z) Radius*SIZE(1,2,3), and the cylinder has the radius Radius*SIZE(4), height Radius*SIZE(5)

[L,M,N,~] = size(r);

Rx=Radius*SIZE(1);
Ry=Radius*SIZE(2);
Rz=Radius*SIZE(3); % for spheroid

Rcy=Radius*SIZE(4);
Hcy=Radius*SIZE(5); % for cylinder

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = ones(L,M,N);

% define spheroid
object1 = @(r)( ((r(:,:,:,1)-Rc(1,1)).^2/Rx^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Ry^2 + (r(:,:,:,3)-Rc(3,1)).^2/Rz^2 < 1) & ((r(:,:,:,3)-Rc(3,1))>0));
pointobject1 = object1(r); %
idx = find(pointobject1(:)); % get indexes of elements
mu_r(idx) = m_r; 
epsilon_r(idx)=e_r(1);
T(idx)=Temp(1);

% define cylinder
object2 = @(r)( ((r(:,:,:,1)-Rc(1,1)).^2/Rcy^2 +  (r(:,:,:,2)-Rc(2,1)).^2/Rcy^2 <= 1) & ((r(:,:,:,3)-Rc(3,1))>-Hcy) ...
              &((r(:,:,:,3)-Rc(3,1))<=0) );
point=object2(r);
idx1=find(point(:));
epsilon_r(idx1) = e_r(2); 
T(idx1) = Temp(2);
end
function [epsilon_r,mu_r,T,idx] = ellipsoid_mix2(r,Rc,Radius,e_r,m_r,Temp,SIZE)
% Function that generates an ellipsoid with semi-axis (major in x,y,z directions): Radius*SIZE(1,2,3). While the part Z>Radiusz-2*Radiusz*SIZE(4) has Temp(1)/e_r(1), and the rest Temp(2)/e_r(2)
    
Radiusx=Radius*SIZE(1);
Radiusy=Radius*SIZE(2);
Radiusz=Radius*SIZE(3); %semi-axis

[L,M,N,~] = size(r);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = ones(L,M,N);

r(:,:,:,1)=r(:,:,:,1)-Rc(1,1);
r(:,:,:,2)=r(:,:,:,2)-Rc(2,1);
r(:,:,:,3)=r(:,:,:,3)-Rc(3,1);
% define upper-part 
object = @(r)( r(:,:,:,1).^2/Radiusx^2 +  r(:,:,:,2).^2/Radiusy^2 + r(:,:,:,3).^2/Radiusz^2 < 1) ;
pointobject = object(r); % ones for domain elements
idx = find(pointobject(:)); % get indexes of elements
mu_r(idx) = m_r; 
epsilon_r(idx) = e_r(1);
T(idx) = Temp(1);

% define lower-part
object1= @(r) (r(:,:,:,3)<=Radiusz-2*Radiusz*SIZE(4));
point1=object1(r);
idx1=find(point1(:)+pointobject(:)==2);
epsilon_r(idx1) = e_r(2); 
T(idx1) = Temp(2);
end
function [epsilon_r,mu_r,T,idx] = hemispheroidshell_Tfromdata(r,Rc,Radius,e_r,m_r,SIZE,dataT,core_diperseT)
%   Function that generates a hemispheroid with shell with temperature distribution from comsol simulation data, and the permittivity of the core is a function of temperature. The shell has the permittivity e_r

%   The size: core has semi-axis SIZE(1,2,3) for x,y,z directions,
%   and the shell has the relative thickness SIZE(4)
        
X=dataT(:,1);Y=dataT(:,2);Z=dataT(:,3);Tem=dataT(:,4);
F=scatteredInterpolant(X,Y,Z,Tem,'linear'); %fit the temperature profile from comsol simulation data

[L,M,N,~] = size(r);

Rx=Radius*SIZE(1);
Ry=Radius*SIZE(2);
Rz=Radius*SIZE(3); %size for spheroid

r(:,:,:,1)=r(:,:,:,1)-Rc(1,1);
r(:,:,:,2)=r(:,:,:,2)-Rc(2,1);
r(:,:,:,3)=r(:,:,:,3)-Rc(3,1);
Xg=r(:,:,:,1);Yg=r(:,:,:,2);Zg=r(:,:,:,3);

epsilon_r = ones(L,M,N);
mu_r = ones(L,M,N);
T = ones(L,M,N);

% define core
core = @(r)( (r(:,:,:,1).^2/Rx^2 +  r(:,:,:,2).^2/Ry^2 + r(:,:,:,3).^2/Rz^2 < 1) & (r(:,:,:,3)>0) );
pointcore = core(r); % ones for domain elements in the shell
idx = find(pointcore(:)); % get indexes of elements
mu_r(idx) = m_r; 

X_core=Xg(idx);Y_core=Yg(idx);Z_core=Zg(idx); %points inside core
T(idx)=F(X_core,Y_core,Z_core);
epsilon_r(idx)=core_diperseT(T(idx));
% define shell part
shell = @(r)( (r(:,:,:,1).^2/Rx^2 +  r(:,:,:,2).^2/Ry^2 + r(:,:,:,3).^2/Rz^2 <= (1+SIZE(4))^2) & ...
              (r(:,:,:,1).^2/Rx^2 +  r(:,:,:,2).^2/Ry^2 + r(:,:,:,3).^2/Rz^2 >= 1) & ...
              (r(:,:,:,3)>0) );
point=shell(r);
idx=find(point(:));
mu_r(idx) = m_r; 
epsilon_r(idx) = e_r(1); 
X_shell=Xg(idx);Y_shell=Yg(idx);Z_shell=Zg(idx);
T(idx)=F(X_shell,Y_shell,Z_shell);
end