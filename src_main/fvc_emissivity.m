function [Power,Ranks] = fvc_emissivity(freq, r, EMT, OPTIONS)
%% Computes the overall thermal radiating energy to the far-field

% voxel volume
Gram = (r(2,1,1,1) - r(1,1,1,1))^3;

% Electromagnetic variables
EM = em_var(freq);

% Get diagonal matrices with the material properties 
[M] = getM2D(EMT.Er,EMT.Tr,EM.ce,Gram,freq);

% Get the operator
fprintf('\n N Operator');
[fN] = getOPERATORS(r,freq,'N');

% Get truncated-SVD
fprintf('\n rSVD -  N Operator');
[U,S] = getUS(fN,M.idx,EM.ce,OPTIONS.SVD_TOL);

% -------------------------------------------------------------------------
% Compute
%         XU = W' * U
% where W = inv(A)
% -------------------------------------------------------------------------

fprintf('\n ------------------------------ \n' )
fprintf('Solve for U of symG \n\n')

XU = zeros(size(U));

% Solve the system for each singular vector
tic
parfor nU = 1:size(U,2)
    fprintf('Solving %3d /%3d \r',nU,size(U,2))
    [XU(:,nU)] = solve_WtransU(U(:,nU),fN,r,EMT,M,OPTIONS);
end
Time_XU = toc;
fprintf('Time_XU      = %dm%ds  \n\n' ,floor(Time_XU/60),int64(mod(Time_XU,60)))

% Square roots of matrices D & S
L_S = sqrtm(S);

L_D = sqrt(diag( M.D .* M.pf ));
L_D = diag(L_D);

% Final results
Power = norm( (L_D' * XU * L_S),'fro')^2;
Ranks = size(S,1);
end