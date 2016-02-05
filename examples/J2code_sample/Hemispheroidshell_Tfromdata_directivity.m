clear
format long
clc
shape = 'Hemispheroidshell_Tfromdata'

micron = 1e-6;
Radius = 0.63*micron;
rx=1;ry=1;rz=1.4; thick=2; %size specification

ObjProperties.Mu           = 1;
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;
ObjProperties.SIZE       = [rx ry rz thick];

load halfellipsincoatT873.txt; %temperature distribution calculated from comsol
ObjProperties.dataT=halfellipsincoatT873;
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
alpha = [1]
freq = alpha * (299792458/(2*pi)) / Radius;

% Discretization
nX = 20;

ObjProperties.Epsilon  = si3n4_permittivity(freq); %permittivity for shell.
ObjProperties.core_diperseT=@(T) GST_permittivity(freq,T) %permittivity for the core.
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



