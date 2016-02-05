function [U,S,V] = rSVD_USV(Fdirect,Fadjoint,m,n,Lmax,tol,blocksize)
%% Computes the truncated (randomized) SVD of an operator

%   INPUT:
%   Fdirect  - handle to direct operator
%   Fadjoint - handle to the adjoint operator
%   m - Left dimension of Fdirect
%   n - right dimension of Fdirect
%   Lmax - maximum number of random excitations (if <1 automatic selection)
%   tol - tolerance for truncation
%   blocksize - block for iterative approaches
%
%   OUTPUT:
%   U - left subspace
%   S - singular values
%   V - right subspace
%
%   ||A - U*S*V'|| / ||A|| < tol
%
% -------------------------------------------------------------------------

if (nargin < 4)
    fprintf(1, '/nERROR: invalid arguments\n');
    U = [];
    S = [];
    V = [];
    return
end

if isempty(Lmax) || (nargin < 5)
    Lmax = min(500,N);
end

if isempty(tol) || (nargin < 6)
    tol = 1e-3;
end

if isempty(blocksize) || (nargin < 7)
    blocksize = 100;
end


% -------------------------------------------------------------------------
%   Random excitations
%   Gm: (mxl)
%   Gn: (nxl)
%
%   Action on random matrices
%   AGn = A*Gn: (mxl) 
%   AGm = A'*Gm: (nxl)
%
%  ON basis
%  [Q,~] = qr(AGn);    
%  [W,~] = qr(AGm);
%
%
%  T = Q'*A*W;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%            Generate the random excitations
% -------------------------------------------------------------------------

fid = 1;
tini = tic;

% if on-the-fly control of SV drop
if (Lmax < 1)
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n RSVD_USV starting with on-the-fly control and tol = %1.1d',tol);
else
    % generate the random excitation matrix, and initialize the fields
    fprintf(fid, '\n\n RSVD_USV starting with %d excitations', Lmax);
end


% -------------------------------------------------------------------------
%            Q and W
% -------------------------------------------------------------------------
fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'ON basis (Q) for left subspace \n');

t2 = tic;
[Q,Nrand_Q] = getQ(Fdirect,m,n,Lmax,tol,blocksize);
timeQR_left = toc(t2);
fprintf(fid,'\n\nelapsed time %g, with %d excitations', timeQR_left,Nrand_Q);

fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'ON basis (W) for right subspace \n');

t2=tic;
[W,Nrand_W] = getQ(Fadjoint,n,m,Lmax,tol,blocksize);
timeQR_right = toc(t2);
fprintf(fid,'\n\nelapsed time %g, with %d excitations', timeQR_right,Nrand_W);

% -------------------------------------------------------------------------
%  T = Q'*A*W;
% -------------------------------------------------------------------------
t1 = tic;

% check minimum columns between Q and W

if size(Q,2) < size(W,2)
    %  compute  A'*Q
    U = zeros(n,size(Q,2));
    for i = 1:size(Q,2)
        
        ve = Q(:,i); 
        
        % adjoint operator
        U(:,i) = Fadjoint(ve);
    end
    % compute T
    T = U'*W;
else
    %  compute  A*W
    U = zeros(m,size(W,2));
    for i = 1:size(W,2)
        
        ve = W(:,i); 
        
        % adjoint operator
        U(:,i) = Fdirect(ve);
    end
    % compute T
    T = Q'*U;
end
timeinc = toc(t1);
fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid,'T computed, elapsed time %g \n', timeinc);

% -------------------------------------------------------------------------
%          [Ur,S,Vr] = svd(T);
% -------------------------------------------------------------------------

t1 = tic;

% -------------------------------------------------------------------------
% apply SVD and truncation with defined tolerance
% -------------------------------------------------------------------------
[U, S, V] = svd(T,'econ');

for kk = 2:length(S)
    if S(kk,kk) < tol*S(1,1)
        Nsv = kk;
        break
    end
end

if kk == length(S)
    Nsv = kk;
end

U = U(:,1:kk);
S = S(1:kk,1:kk);
V = V(:,1:kk);
% -------------------------------------------------------------------------
% Final output
% -------------------------------------------------------------------------
U = Q*U; 
V = W*V;

timeinc = toc(t1);
fprintf(fid,'\n SVD(T) done, time %g', timeinc);

fprintf(fid, '\n\n-------------------------- \n');
fprintf(fid, 'Random SVD operation done\n');
fprintf(fid, ' elapsed time %g seconds\n',toc(tini));
fprintf(fid, ' #Singular values  %d\n',Nsv);

end