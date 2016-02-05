function [xJ] = solve_WEinc(Ein,fN,R,EMT,M,freq,OPTIONS)
%% Solve the JVIE linear system: A*x = ce*chi*Einc => x = W*(ce*chi*Einc)

fid = 1; % for fprintf

% EM variables
EM =  em_var(freq);

% voxel volume
dx = R(2,1,1,1) - R(1,1,1,1);
dV = dx^3;

% -------------------------------------------------------------------------
% Excitation pre-processing
% -------------------------------------------------------------------------

% translate the input in the dimensions of the scatterer
% set the multiplier to convert fields to suitable input
sE = zeros(length(EMT.Er(:)),3);
for ii = 1:3
    vE = (EM.ce*(EMT.Er-1)) .* Ein(:,:,:,ii);
    sE(:,ii) = vE(:);
end
Vexc = sE(:);
clear vE sE

% right hand side
b = Vexc(M.idx);
clear Vexc

% define the main function
fA   = @(J)mv_AN(J, fN, EMT.Er, EMT.Er-1, dV, 1, M.idx, 'notransp');

% Preconditioner
switch OPTIONS.PRECOND
    case 0
        fP = [];
    case 1
        fP =  M.Mer;
    case 2
        fP =  @(J)mv_AN(J, fN, 1./EMT.Er, (1./EMT.Er)-1, dV, 1, M.idx, 'notransp');
end
       
% -------------------------------------------------------------------------
% Select the iterative method and compute the inverse
% -------------------------------------------------------------------------
   
if (OPTIONS.ITSOLVER == 1) %  BICGSTAB
     
    tsolve = tic;
    [xJ,flag_,~,~,resvec_] = bicgstab(@(J)fA(J), b, OPTIONS.TOL, OPTIONS.INNER_IT*OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = b - fA(xJ);
        fprintf(fid,'\n MATLAB BUILT IN BICGSTAB FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(b));
    end
           
else %  GMRES
    
    tsolve = tic;
    [xJ,flag_,~,~,resvec_] = gmres(@(J)fA(J), b, OPTIONS.INNER_IT, OPTIONS.TOL, OPTIONS.OUTER_IT,fP);
    solve_time = toc(tsolve);
    
    if (OPTIONS.VERBOSE == 1)
        r_ = b - fA(xJ);
        fprintf(fid,'\n MATLAB BUILT IN GMRES FUNCTION, NO PRECONDITIONER\n');
        fprintf(fid,' Time  = %g [sec], flag %d, %d iterations, residue %g, relative residue %g \n' ,solve_time,flag_,length(resvec_),norm(r_),norm(r_)/norm(b));
    end
           
end

end



