function [U,S] = getUS(fN,idx,ce,eig_tol)
%% Computes the truncated (randomized) SVD of (symmetric) matrix symG

% operator symG
Fdirect = @(x)mv_symG(x, fN, idx, ce);

% input for rSVD
Mrsvd = length(idx);
Nrsvd = length(idx);
Lmax = 1e-1;
tol  = eig_tol;
blocksize = 10;

% main function
[U,S] = rSVD_US(Fdirect,Fdirect,Mrsvd,Nrsvd,Lmax,tol,blocksize);

end
    