function [xU] = solve_WtransU_2b(U,fN_11,fN_12,EMT_1, EMT_2, M1, M2, Gram, OPTIONS, mode)
%% Solve the JVIE-based linear system for 2 bodies: ctranspose(A)*x = U => x = ctranspose(W)*U

fid = 1;

% select "mode" and right hand side
sz1 = length(M1.idx);
sz2 = length(M2.idx);

if strcmp(mode,'11')
    
    assert(size(U,1) == sz1, 'size mismatch')
    
    Ub = [U;zeros(sz2,1)];
    
elseif strcmp(mode,'12')
    
    assert(size(U,1) == sz2, 'size mismatch')
    
    Ub = [zeros(sz1,1);U];
    
elseif strcmp(mode,'21')
    
    assert(size(U,1) == sz1, 'size mismatch')
    
    Ub = [U;zeros(sz2,1)];
    
elseif strcmp(mode,'22')
    
    assert(size(U,1) == sz2, 'size mismatch')
    
    Ub = [zeros(sz1,1);U];
end

% Define the main function
fA = @(J)mv_AN_2b(J, fN_11, fN_12, EMT_1, EMT_2, M1, M2, Gram, 'transp');

% Preconditioner
switch OPTIONS.PRECOND
    case 0
        fP = [];
    case 1
        fP =  blkdiag(M1.Mer,M2.Mer);
    case 2
        EMT_1.Er = 1./EMT_1.Er;
        EMT_2.Er = 1./EMT_2.Er;
        %
        fP = @(J)mv_AN_2b(J, fN_11, fN_12, EMT_1, EMT_2, M1, M2, Gram, 'transp');
end
        
% -------------------------------------------------------------------------
% Select the iterative method and compute the inverse
% -------------------------------------------------------------------------
   
if (OPTIONS.ITSOLVER == 1) %  bICGSTAb

    tsolve = tic;
    [vsol_,flag_,~,~,resvec_] = bicgstab(@(J)fA(J), Ub, OPTIONS.TOL, OPTIONS.INNER_IT*OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = Ub - fA(vsol_);
        fprintf(fid,'\n MATLAb bUILT IN bICGSTAb FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(U));
    end
     
else %  GMRES

    tsolve = tic;
    [vsol_,flag_,~,~,resvec_] = gmres(@(J)fA(J), Ub, OPTIONS.INNER_IT, OPTIONS.TOL, OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = Ub - fA(vsol_);
        fprintf(fid,'\n MATLAB BUILT IN GMRES FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(U));
    end
  
end

% Output
if strcmp(mode,'11')
    
    xU = vsol_(1:sz1,1);
    
elseif strcmp(mode,'12')
    
    xU = vsol_(1:sz1,1);
    
elseif strcmp(mode,'21')
    
    xU = vsol_(sz1+1:end,1);
    
elseif strcmp(mode,'22')
    
    xU = vsol_(sz1+1:end,1);
end

end


