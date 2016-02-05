function LDOSU = fvc_LDOSdirectivity(DIRECTIONS,freq,r,EMT,OPTIONS)
%% Compute the thermal radiation from each volume for each direction in far-field

% voxel volume
Gram = (r(2,1,1,1) - r(1,1,1,1))^3;

% Electromagnetic variables
EM = em_var(freq);

% Get the operator
[fN] = getOPERATORS(r,freq,'N');

% Get diagonal matrices with the material properties 
[M] = getM2D(EMT.Er,EMT.Tr,EM.ce,Gram,freq);

% Compute radiation intensity for each direction
LDOSU = getULDOS(fN,DIRECTIONS,freq,r,EMT,M,OPTIONS);

end
%% Local functions
function [UL] = getULDOS(fN,DIRECTIONS,f,r,EMT,M,OPTIONS)

n_t = size(DIRECTIONS.theta,1);
n_p = size(DIRECTIONS.phi,1);

assert(n_t == n_p, 'Theta and Phi should have the same dimension')

% -------------------------------------------------------------------------
% Get directivity intensity
% -------------------------------------------------------------------------

% First compute matrices H1, H2

[Htp1,Htp2] = getH(DIRECTIONS,f,r,EMT.Er,M.idx);

% Get the conjugate transpose

Htp1 = ctranspose(Htp1);
Htp2 = ctranspose(Htp2);

[UL] = u_H(Htp1,Htp2,f,r,fN,EMT,M,OPTIONS);

end
function [U] = u_H(H1,H2,freq,r,fN,EMT,M,OPTIONS)

% EM variables
EM = em_var(freq); 

% pre-factor
U_pre = EM.ko^4  / ((4*pi)^2 * EM.eta * abs(EM.ce)^2);  % 1/2 included in definition of M.D!

% -------------------------------------------------------------------------
% Correlation matrix
% -------------------------------------------------------------------------
L_D = sqrt(diag( M.D .* M.pf ));
L_D = diag(L_D);

% number of directions
nL = size(H1,2);

% -------------------------------------------------------------------------
% Trace computation
%    P = Tr[(L_D * W' * H' * P23) * (L_D * W' * H' * P23)']
%      = Tr[Q * Q']
%      = |Q|^2
% -------------------------------------------------------------------------

H = [H1 H2]; % this is H' !!
U_fro = zeros(length(M.idx),size(H,2));

% for -> system solution for all vectors (directions * 2)
check_matlabpool;

parfor ctL = 1: 2*nL
    
    % -------------------------------------------------------------------------
    % Compute xU = W'*U -> xU = [A^(-1)]' * U
    % -------------------------------------------------------------------------
    
    xH = solve_WtransU(H(:,ctL),fN,r,EMT,M,OPTIONS);
    
    xH = L_D * xH;
    U_fro(:,ctL) = abs(xH).^2;
end
U = U_fro(:,1:nL) + U_fro(:,nL+1:2*nL);

% Final output
U = U_pre * U;

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

end 
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