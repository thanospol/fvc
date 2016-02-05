clear
clc
format long
% -------------------------------------------------------------------------
% EXAMPLE - GEOMETRY
% -------------------------------------------------------------------------
% shape = 'Sphere'
%       = 'Cube'
%       = 'Cylinder'

shape = 'Cube';

micron = 1e-6;
Radius = micron;
Diam = 2*Radius;

ObjProperties.Epsilon      = 10-1i*1;
ObjProperties.Mu           = 1;
ObjProperties.Temp_profile = 300;
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;

% -------------------------------------------------------------------------
% OPTIONS - for the iterative solver and the SVD_tol
% -------------------------------------------------------------------------

OPTIONS.ITSOLVER    = 1; % ITSOLVER (1) for BICGSTAB, GMRES otherwise
OPTIONS.TOL         = 1e-3; % TOL tolerance for solver
OPTIONS.OUTER_IT    = 50; % OUTER_IT outer iterations for GMRES (for BICGSTAB the overall number is INNER_IT*OUTER_IT)
OPTIONS.INNER_IT    = 100; % INNER_IT inner iterations for GMRES
OPTIONS.VERBOSE     = 0; % VERBOSE (1) if wanna see print info, no printing as default
OPTIONS.PRECOND     = 0; % PRECOND: (0) no preconditioner
                         %          (1) left preconditioner for highly inhomogeneous objects
                         %          (2) left preconditioner for high contrast    
OPTIONS.SVD_TOL     = 1e-3; % SVD_TOL tolerance for truncated SVD


% -------------------------------------------------------------------------
% INPUT 
% -------------------------------------------------------------------------

% Discretization
nX = 11;

[r,EMT] = getGeometry_1obj(nX,ObjProperties);

% Frequency range
alpha = [1.0];
freq = alpha * (299792458/(2*pi)) / ObjProperties.Radius;

% -------------------------------------------------------------------------
% RUN fvc_emissivity
% -------------------------------------------------------------------------

[Power,Ranks] = fvc_emissivity(freq, r, EMT, OPTIONS);

Power
