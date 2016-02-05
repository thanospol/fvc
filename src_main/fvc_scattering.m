function [U,xJ,Etot] = fvc_scattering(freq,r,EMT,OPTIONS,DIRECTIONS)
%% Computes the far and local field intensity after scattering.
% the far-field intensity will be calculated if DIRECTIONS is given.

% voxel volume
Gram = (r(2,1,1,1) - r(1,1,1,1))^3;

% Electromagnetic variables
EM = em_var(freq);

% Get the operators
[fN] = getOPERATORS(r,freq,'N');
%[fK] = getOPERATORS(r,freq,'K'); % needed for H-field

% Get diagonal matrices with the material properties 
[M] = getM2D(EMT.Er,EMT.Tr,EM.ce,Gram,freq);

% -------------------------------------------------------------------------
% Solve the scattering problem
% -------------------------------------------------------------------------

% Define excitation - a plane wave
polarization = 'x';
k_vector = [0,0,EM.ko]; % direction
[Einc, ~] = planewave(r,k_vector,EM.omega_mu,polarization);

% -------------------------------------------------------------------------
% FFT VIE solver
% -------------------------------------------------------------------------
[xJ] = solve_WEinc(Einc,fN,r,EMT,M,freq,OPTIONS);

% -------------------------------------------------------------------------
% EH-fields from J-currents
% -------------------------------------------------------------------------
[Etot] = e_field_Nop_comp(xJ, fN, Gram, EM.omega, EM.eo, Einc, M.idx);

% -------------------------------------------------------------------------
% Get (Scattered) radiation intensity
% -------------------------------------------------------------------------

if nargin==4
    U=[];
elseif nargin==5
    U = getUsca(xJ,DIRECTIONS,freq,r,EMT.Er,M.idx);
end

end
%% Local functions
function [UL_sca] = getUsca(xJ,DIRECTIONS,freq,r,Er,idx)
%%
n_t = size(DIRECTIONS.theta,1);
n_p = size(DIRECTIONS.phi,1);

assert(n_t == n_p, 'Theta and Phi should have the same dimension')

nL = n_t;

% EM variables
EM = em_var(freq);

% pre-factor
U_pre = EM.ko^4  / (2 * (4*pi)^2 * EM.eta * abs(EM.ce)^2);  
% output vector
UL_sca = zeros(nL,1);

% -------------------------------------------------------------------------
% Get directivity intensity
% -------------------------------------------------------------------------

% First compute matrices H1, H2

[H1,H2] = getH(DIRECTIONS,freq,r,Er,idx);

% Compute H*x 
xH_1 = H1 * xJ;
xH_2 = H2 * xJ;

for ctL = 1: nL

        % -------------------------------------------------------------------------
        % Compute xU = W'*U -> xU = [A^(-1)]' * U
        % -------------------------------------------------------------------------

        xH = [xH_1(ctL,1);xH_2(ctL,1)];
        
        UL_sca(ctL,1) = xH' * xH;

end

% Final output
UL_sca =  U_pre * UL_sca;
end
function [H1,H2] = getH(DIRECTIONS,freq,r,Er,idx)
%% evaluate matrix H
n_t = size(DIRECTIONS.theta,1);
n_p = size(DIRECTIONS.phi,1);

assert(n_t == n_p, 'Theta and Phi should have the same dimension')

nL = n_t;

% dx
dx = r(2,1,1,1)-r(1,1,1,1);

% EM variables
EM = em_var(freq);

% pre-factor

P23 = [0 0 0 
       0 1 0
       0 0 1];
   
% Pre-allocation 
F  = zeros(nL,1);
H1 = zeros(nL,length(idx));
H2 = zeros(nL,length(idx));

for ctL = 1: nL
    
    theta = DIRECTIONS.theta(ctL);
    phi   = DIRECTIONS.phi(ctL);
    
    Q_C2S = [sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta);
             cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta);
            -sin(phi),            cos(phi),             0.0 ];
    
    Q = P23 * Q_C2S;
    
    % Evaluate scalar f
    
    Ftp = getFtp(EM.ko,theta,phi,dx);
    F(ctL) = Ftp;
    
    % Form matrix H
    
    H = geth_(EM.ko,theta,phi,r,Er);
    
    H = Q * H;
    
    % We only need the 2 non-zero columns
    
    H1(ctL,:) = H(2,:);
    H2(ctL,:) = H(3,:);
    
    %
    H1(ctL,:) = F(ctL,1) * H1(ctL,:);
    H2(ctL,:) = F(ctL,1) * H2(ctL,:);
end

end % function [H1,H2] = GetH(DIRECTIONS,freq,r,Er,idx)
function Ftp = getFtp(ko,inclination,azimuth,Dx)

D2 = Dx/2;
%

x = sin(inclination)*cos(azimuth);
if abs(x)>eps
    Ftp_1 = 2 * sin(ko*x*D2) / (ko*x);
else
    Ftp_1 = 2*D2;
end

x = sin(inclination)*sin(azimuth);
if abs(x)>eps
    Ftp_2 = 2 * sin(ko*x*D2) / (ko*x);
else
    Ftp_2 = 2*D2;
end

x = cos(inclination);
if abs(x)>eps
    Ftp_3 = 2 * sin(ko*x*D2) / (ko*x);
else
    Ftp_3 = 2*D2;
end
% Output
Ftp = Ftp_1 * Ftp_2 * Ftp_3;
end
function hmatrix = geth_(ko,inclination,azimuth,r,Er)

% univec_x = [ sin(inclination)*cos(phi); sin(inclination)*sin(phi); cos(inclination) ];

x_x = sin(inclination)*cos(azimuth);
x_y = sin(inclination)*sin(azimuth);
x_z = cos(inclination);

rx = r(:,:,:,1);
ry = r(:,:,:,2);
rz = r(:,:,:,3);

R = rx*x_x + ry*x_y + rz*x_z;

% BE CAREFUL: here it is exp(+1i*ko*R)!!!
h = exp(+1i*ko*R);

% Get the values where there is material
hh = h((Er-1.0)~=0);

hh = hh(:);

hh = transpose(hh);

% Form the final matrix
hmatrix = blkdiag(hh,hh,hh);
end