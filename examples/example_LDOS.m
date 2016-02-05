clear
format long
clc
% -------------------------------------------------------------------------
% EXAMPLE - GEOMETRY
% -------------------------------------------------------------------------
% shape = 'Sphere'
%       = 'Cube'
%       = 'Cylinder'

shape = 'Sphere';

micron = 1e-6;
Radius = micron;
Diam = 2*Radius;

ObjProperties.Epsilon      = 10-1i*1;
ObjProperties.Mu           = 1;
ObjProperties.Temp_profile = 300;
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;

optionplot=1; %1 for plotting LDOS
% -------------------------------------------------------------------------
% OPTIONS - for the iterative solver 
% -------------------------------------------------------------------------
OPTIONS.ITSOLVER    = 2; % ITSOLVER (1) for BICGSTAB, GMRES otherwise
OPTIONS.TOL         = 1e-3; % TOL tolerance for solver
OPTIONS.OUTER_IT    = 50; % OUTER_IT outer iterations for GMRES (for BICGSTAB the overall number is INNER_IT*OUTER_IT)
OPTIONS.INNER_IT    = 100; % INNER_IT inner iterations for GMRES
OPTIONS.VERBOSE     = 1; % VERBOSE (1) if wanna see print info, no printing as default
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
alpha = [1];
freq = alpha * (299792458/(2*pi)) / ObjProperties.Radius;

% -------------------------------------------------------------------------
% Lebedev Quadrature rule
% -------------------------------------------------------------------------
   
Lebedev_degree = 38;

[w, phiL, thetaL] = lebedev(Lebedev_degree);

DIRECTIONS.theta = thetaL;
DIRECTIONS.phi   = phiL;

% -------------------------------------------------------------------------
% Get radiation intensity from each volume
% -------------------------------------------------------------------------

LDOSU = fvc_LDOSdirectivity(DIRECTIONS,freq,r,EMT,OPTIONS);

LDOS = sum(LDOSU*w,2);

if optionplot
    plot_LDOS(r,EMT.Er,Radius,LDOS);
end


