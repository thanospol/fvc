function [Power,Ranks] = fvc_pt2bodies(freq, r_1, r_2, r_offset, EMT_1, EMT_2, OPTIONS)
%% INITIALIZATION

% voxel volume
Gram = (r_1(2,1,1,1) - r_1(1,1,1,1))^3;

% Electromagnetic variables
EM = em_var(freq);

% Get matrices M
M1 = getM2D(EMT_1.Er,EMT_1.Tr,EM.ce,Gram, freq);
M2 = getM2D(EMT_2.Er,EMT_2.Tr,EM.ce,Gram, freq);

% -------------------------------------------------------------------------
% Get the operators
% -------------------------------------------------------------------------

% self block
fprintf('\n N11 Operator');
[fN_11] = getOPERATORS(r_1,freq,'N');

% coupling block
fprintf('N12 Operator');
[fN_12] = getOPERATORS_C(r_1,-r_offset,freq,'N');

% -------------------------------------------------------------------------
% Get truncated-SVD
% -------------------------------------------------------------------------

fprintf('\n rSVD -  N11 Operator');
[U11,S11] = getUS(fN_11,M1.idx,EM.ce,OPTIONS.SVD_TOL);

fprintf('\n rSVD -  N12 Operator');
[U12,S12,V12] = getUSV(fN_12,M1.idx,EM.ce,OPTIONS.SVD_TOL);

% -------------------------------------------------------------------------
% [~;XU11] = W' * [U11;0]
% -------------------------------------------------------------------------
fprintf('\n ------------------------------ \n' )
fprintf('Solve for U11 of symG \n\n')

XU11 = zeros(size(U11));

% Solve the system for each singular vector
tic
parfor nU11 = 1:size(U11,2)
    fprintf('Solving %3d /%3d \r',nU11,size(U11,2))
    [XU11(:,nU11)] = solve_WtransU_2b(U11(:,nU11),fN_11,fN_12,EMT_1, EMT_2, M1, M2, Gram, OPTIONS,'21')
end
Time_XU11 = toc;
fprintf('Time_XU11      = %dm%ds  \n\n' ,floor(Time_XU11/60),int64(mod(Time_XU11,60)))

% -------------------------------------------------------------------------
% [~;XU12] = W' * [U12;0]
% -------------------------------------------------------------------------
fprintf('\n ------------------------------ \n' )
fprintf('Solve for U12 of G12 \n\n')

XU12 = zeros(size(U12));

% Solve the system for each singular vector
tic
parfor nU12 = 1:size(U12,2)
    fprintf('Solving %3d /%3d \r',nU12,size(U12,2))
    [XU12(:,nU12)] = solve_WtransU_2b(U12(:,nU12),fN_11,fN_12,EMT_1, EMT_2, M1, M2, Gram, OPTIONS, '21')
end
Time_XU12 = toc;
fprintf('Time_XU12      = %dm%ds  \n\n' ,floor(Time_XU12/60),int64(mod(Time_XU12,60)))

% -------------------------------------------------------------------------
% [~;XV12] = W' * [0;V12]
% -------------------------------------------------------------------------
fprintf('\n ------------------------------ \n' )
fprintf('Solve for V12 of G12 \n\n')

XV12 = zeros(size(V12));

% Solve the system for each singular vector
tic
parfor nV12 = 1:size(V12,2)
    fprintf('Solving %3d /%3d \r',nV12,size(V12,2))
    [XV12(:,nV12)] = solve_WtransU_2b(V12(:,nV12),fN_11,fN_12,EMT_1, EMT_2, M1, M2, Gram, OPTIONS, '22')
end
Time_XV12 = toc;
fprintf('Time_XV12      = %dm%ds  \n\n' ,floor(Time_XV12/60),int64(mod(Time_XV12,60)))

% -------------------------------------------------------------------------
% Square roots of matrices D & S
% -------------------------------------------------------------------------
L_S11 = sqrtm(S11);
L_S12 = sqrtm(S12);

L_D2 = chol(M2.D);

% -------------------------------------------------------------------------
% Traces
% -------------------------------------------------------------------------

% "self-term", Frobenius norm
P21_11 = norm( (L_D2' * XU11 * L_S11),'fro')^2;

% "coupling" term
XU12 = L_D2' * XU12 * L_S12;
XU12 = XU12(:);
XU12 = transpose(XU12);

XV12 = L_D2' * XV12 * L_S12;
XV12 = XV12(:);
XV12 = conj(XV12);

P21_12 = XU12 * XV12;
P21_12 = real(P21_12);

% power transfer
P21 = P21_12 - P21_11;

% -------------------------------------------------------------------------
% Collect results
% -------------------------------------------------------------------------
Power.self     = P21_11;
Power.coupling = P21_12;
Power.transfer = P21;

Ranks.N11 = size(S11,1);
Ranks.N12 = size(S12,1);
end