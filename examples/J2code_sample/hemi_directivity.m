clear
format long
clc
shape = 'Hemispheroidshell';

micron = 1e-6;
Radius = micron;
rx=1;ry=1;rz=.4;thick=.4; %size specification

ObjProperties.Epsilon      = [30-5i -100-80i];
ObjProperties.Mu           = 1;
ObjProperties.Temp_profile = [1000 0];
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;
ObjProperties.SIZE       = [rx ry rz thick];

% -------------------------------------------------------------------------
% OPTIONS - for the iterative solver 
% -------------------------------------------------------------------------
OPTIONS.ITSOLVER    = 2; % ITSOLVER (1) for BICGSTAB, GMRES otherwise
OPTIONS.TOL         = 1e-4; % TOL tolerance for solver
OPTIONS.OUTER_IT    = 50; % OUTER_IT outer iterations for GMRES (for BICGSTAB the overall number is INNER_IT*OUTER_IT)
OPTIONS.INNER_IT    = 100; % INNER_IT inner iterations for GMRES
OPTIONS.VERBOSE     = 1; % VERBOSE (1) if wanna see print info, no printing as default
OPTIONS.PRECOND     = 2; % PRECOND: (0) no preconditioner
                         %          (1) left preconditioner for highly inhomogeneous objects
                         %          (2) left preconditioner for high contrast    
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------

% Frequency range
alpha = [.8]
freq = alpha * (299792458/(2*pi)) / Radius;

% Discretization
nX = 20;
[r,EMT] = getGeometry_1obj(nX,ObjProperties);
% -------------------------------------------------------------------------
% Lebedev Quadrature rule
% -------------------------------------------------------------------------
   
Lebedev_degree = 14;

[w, phiL, thetaL] = lebedev(Lebedev_degree);

DIRECTIONS.theta = thetaL;
DIRECTIONS.phi   = phiL;

% -------------------------------------------------------------------------
% Get radiation intensity
% -------------------------------------------------------------------------

U = fvc_directivity(DIRECTIONS,freq,r,EMT,OPTIONS);



