clear
format long
clc
% -------------------------------------------------------------------------
% EXAMPLE - GEOMETRY
% -------------------------------------------------------------------------
% shape = 'Sphere'
%       = 'Cube'
%       = 'Cylinder'

shape = 'Cubelinearz';

micron = 1e-6;
Radius = micron;
Diam = 2*Radius;

ObjProperties.Epsilon      = [2-1i*1 12-1i*1];
ObjProperties.Mu           = 1;
ObjProperties.Temp_profile = [1000 1000];
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;

% -------------------------------------------------------------------------
% OPTIONS - for the iterative solver 
% -------------------------------------------------------------------------
OPTIONS.ITSOLVER    = 2; % ITSOLVER (1) for BICGSTAB, GMRES otherwise
OPTIONS.TOL         = 1e-3; % TOL tolerance for solver
OPTIONS.OUTER_IT    = 50; % OUTER_IT outer iterations for GMRES (for BICGSTAB the overall number is INNER_IT*OUTER_IT)
OPTIONS.INNER_IT    = 100; % INNER_IT inner iterations for GMRES
OPTIONS.VERBOSE     = 0; % VERBOSE (1) if wanna see print info, no printing as default
OPTIONS.PRECOND     = 2; % PRECOND: (0) no preconditioner
                         %          (1) left preconditioner for highly inhomogeneous objects
                         %          (2) left preconditioner for high contrast    
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------

% Discretization
nX = 30;

[r,EMT] = getGeometry_1obj(nX,ObjProperties);

% Frequency range
alpha_array = .1:.1:2;


% -------------------------------------------------------------------------
% Lebedev Quadrature rule
% -------------------------------------------------------------------------
   
Lebedev_degree = 38;

[w, phiL, thetaL] = lebedev(Lebedev_degree);

DIRECTIONS.theta = thetaL;
DIRECTIONS.phi   = phiL;

% -------------------------------------------------------------------------
% Get radiation intensity
% -------------------------------------------------------------------------

ctrl=0;
Prad=zeros(size(alpha_array));
for alpha=alpha_array
    ctrl=ctrl+1;
    freq = alpha * (299792458/(2*pi)) / ObjProperties.Radius;
    U = fvc_directivity(DIRECTIONS,freq,r,EMT,OPTIONS);
    Prad(ctrl) = sum(w.*U);
end

% -------------------------------------------------------------------------
% Get emissivity
% -------------------------------------------------------------------------
C=4.79924335e-11; %C=h/K_b
h=6.62606957e-34;

Freq = alpha_array * (299792458/(2*pi)) / ObjProperties.Radius;
A=6*(2*ObjProperties.Radius)^2; %surface area for this prolate spheroid.
plank = h*Freq./(exp(C*Freq/ObjProperties.Temp_profile(1))-1); %Planck function
black=Freq.^2.*plank/(299792458)^2*A; %Blackbody radiation
emissivity=Prad./black;
