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

% distance between the 2 bodies

d_factor = 2.0;
size_d = length(d_factor);

r_offset = [Diam + d_factor*Diam , 0, 0];

% -------------------------------------------------------------------------
% OPTIONS - for the iterative solver and the SVD_tol
% -------------------------------------------------------------------------

OPTIONS.ITSOLVER    = 2; % ITSOLVER (1) for BICGSTAB, GMRES otherwise
OPTIONS.TOL         = 1e-2; % TOL tolerance for solver
OPTIONS.OUTER_IT    = 50; % OUTER_IT outer iterations for GMRES (for BICGSTAB the overall number is INNER_IT*OUTER_IT)
OPTIONS.INNER_IT    = 100; % INNER_IT inner iterations for GMRES
OPTIONS.VERBOSE     = 0; % VERBOSE (1) if wanna see print info, no printing as default
OPTIONS.PRECOND     = 0; % PRECOND: (0) no preconditioner
                         %          (1) left preconditioner for highly inhomogeneous objects
                         %          (2) left preconditioner for high contrast 
OPTIONS.SVD_TOL     = 1e-2; % SVD_TOL tolerance for truncated SVD

% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------

% Discretization
nX = 11;

[r_1,r_2,EMT_1,EMT_2] = getGeometry_2obj(nX,shape,Radius,r_offset);

% Frequency range
alpha = 1.0;
freq = alpha * (299792458/(2*pi)) / Radius;


[P,R] = fvc_pt2bodies(freq, r_1, r_2, r_offset, EMT_1, EMT_2, OPTIONS);

