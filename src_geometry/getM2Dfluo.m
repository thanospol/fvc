function [M] = getM2Dfluo(er,t,ce,Gram)
%% Get diagonal matrices with the material properties for fluorescence problem.
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
%emitter distribution
pf1 = [t(:);t(:);t(:)];
pf2 = pf1(idx_obj);
M.pf    = spdiags( pf2, 0, Ndof_obj, Ndof_obj );
