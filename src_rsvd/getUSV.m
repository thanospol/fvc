function [U,S,V] = getUSV(fN,idx,ce,eig_tol)
%% Computes the truncated (randomized) SVD of matrix A

% operator N 
Fdirect  = @(x)mv_AN(x, fN, 0, -1, 0, ce, idx, 'notransp');
Fadjoint = @(x)mv_AN(x, fN, 0, -1, 0, ce, idx, 'transp');

% input for rSVD
Mrsvd = length(idx);
Nrsvd = length(idx);
Lmax = 1e-1;
tol  = eig_tol;
blocksize = 50;

% main function
[U,S,V] = rSVD_USV(Fdirect,Fadjoint,Mrsvd,Nrsvd,Lmax,tol,blocksize);

end