function [xJ] = solve_WtransU(b,fN,R,EMT,M,OPTIONS)
%% Solve the JVIE-based linear system: ctranspose(A)*x = b => x = ctranspose(W)*b

fid = 1; % for fprintf

% voxel volume
dx = R(2,1,1,1) - R(1,1,1,1);
dV = dx^3;

% Define the main function
fA   = @(J)mv_AN(J, fN, EMT.Er, EMT.Er-1, dV, 1, M.idx, 'transp');

% Preconditioner
switch OPTIONS.PRECOND
    case 0
        fP = [];
    case 1
        fP =  M.Mer;
    case 2
        fP =  @(J)mv_AN(J, fN, 1./EMT.Er, (1./EMT.Er)-1, dV, 1, M.idx, 'transp');
end
       
% -------------------------------------------------------------------------
% Select the iterative method and compute the inverse
% -------------------------------------------------------------------------
   
if (OPTIONS.ITSOLVER == 1) %  BICGSTAB

    tsolve = tic;
    [xJ,flag_,~,~,resvec_] = bicgstab(fA, b, OPTIONS.TOL, OPTIONS.INNER_IT*OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = b - fA(xJ);
        fprintf(fid,'\n MATLAB BUILT IN BICGSTAB FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(b));
    end
    
else %  GMRES

    tsolve = tic;
    [xJ,flag_,~,~,resvec_] = gmres(fA, b, OPTIONS.INNER_IT, OPTIONS.TOL, OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = b - fA(xJ);
        fprintf(fid,'\n MATLAB BUILT IN GMRES FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(b));
    end
       
end

end


