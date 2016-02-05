function [M] = getM2D(er,t,ce,Gram,freq)
%% Get diagonal matrices with the material properties.

d_chiEx  = er(:) - 1.0;
d_chiE = [d_chiEx;d_chiEx;d_chiEx];

idx_obj = find(d_chiE~=0);
Ndof_obj = length(idx_obj);

% chiE = zeros(Ndof_obj);
chiE = d_chiE(idx_obj);
% size(chiE)

chiE_inv = zeros(size(d_chiE));
% size(chiE_inv)
chiE_inv(idx_obj) = 1./d_chiE(idx_obj);

% Nt = 3*length(d_chiEx(:));

M.MchiE    = spdiags( chiE, 0, Ndof_obj, Ndof_obj );
M.invMchiE = spdiags( chiE_inv, 0, Ndof_obj, Ndof_obj );

M.Id = speye(Ndof_obj,Ndof_obj);
M.Mer = M.MchiE + M.Id;
% idx = find(MchiE);
M.idx = idx_obj;
M.dof = length(idx_obj);

% Planck function
C=4.79924335e-11; %C=h/k_b
h=6.62606957e-34;
pf=h*freq./(exp(C*freq./t)-1);

pf1 = [pf(:);pf(:);pf(:)];
pf2 = pf1(idx_obj);
M.pf    = spdiags( pf2, 0, Ndof_obj, Ndof_obj );

% D (correalation matrix)
M.D = (4/pi) * ce * M.MchiE;
M.D = (M.D + M.D')/2;
M.D = Gram * M.D;
M.D = (1/2) * M.D; % pre-factor 1/2 absorbed here