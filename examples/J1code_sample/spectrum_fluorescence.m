clear;
format long
clc
% -------------------------------------------------------------------------
% EXAMPLE - GEOMETRY
% -------------------------------------------------------------------------
% shape = 'Sphere'
%       = 'Cube'
%       = 'Cylinder'

shape = 'Sphere';
Nphoton=1; %1-photon process

micron = 1e-6;
Radius = micron;
Diam = 2*Radius;

ObjProperties.Epsilon      = 12-1i;
ObjProperties.Mu           = 1;
ObjProperties.Temp_profile = 1; %dipole homogeneous distribution
ObjProperties.shape        = shape;
ObjProperties.Radius       = Radius;

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
nX = 40;

[r,EMT] = getGeometry_1obj(nX,ObjProperties);

% Frequency range
alphaex = 1.58;freqex = alphaex * (299792458/(2*pi)) / ObjProperties.Radius;
alphafluo = 0.1:.1:1.57;
% -------------------------------------------------------------------------
% Lebedev Quadrature rule
% -------------------------------------------------------------------------
   
Lebedev_degree = 38;

[w, phiL, thetaL] = lebedev(Lebedev_degree);
        
DIRECTIONS.theta = thetaL;
DIRECTIONS.phi   = phiL;


% -------------------------------------------------------------------------
% Get (scattered) radiation intensity
% -------------------------------------------------------------------------
[~,xJ,Etot] = fvc_scattering(freqex,r,EMT,OPTIONS);

ctrl=0;
Prad=zeros(size(alphafluo));
for alpha=alphafluo
    ctrl=ctrl+1;
    freq = alpha * (299792458/(2*pi)) / ObjProperties.Radius;
    U = fvc_fluorescence(DIRECTIONS,freq,r,EMT,OPTIONS,Nphoton,Etot);
    Prad(ctrl) = sum(w.*U);
end

      

                


