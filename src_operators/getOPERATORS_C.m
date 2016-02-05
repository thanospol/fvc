function [op_out] = getOPERATORS_C(r,r_offset,f,op,type)
%% function that retrieves the coupling operators N12 and K12 for JM-VIE
%
%% INPUT
% r         3-D Domain
% f         working frequency
%
% op        flag for choosing the operator
%                N: only operator N
%                K: only operator K
%
%% OPTIONAL INPUT
% type      the type for the chosen operator
%                empty: FFT of the Circulant (DEFAULT)
%                T: Toeplitz
%                C: Circulant
% singular  the method for the singular integrals in operator N
%                empty: DIRECTFN (DEFAULT)
%                DEMCEM
%
%% OUTPUT
% op_out    discrete operator
%
% -------------------------------------------------------------------------
%
%   A. G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 5 || isempty(type))
   type = 'FFTofC';
end


% -------------------------------------------------------------------------
% Domain & Frequency
% -------------------------------------------------------------------------
dx = r(2,1,1,1) - r(1,1,1,1); % voxel edge length
ko = 2*pi / 299792458 * f;    % wavenumber, ko = 2*pi/lambda = 2*pi/(co/f);
[Lr,Mr,Nr,~] = size(r);

% -------------------------------------------------------------------------
% Parallelization
% -------------------------------------------------------------------------
if verLessThan('matlab', '8.3.0')
    isOpen = (matlabpool('size') > 0); %#ok<*DPOOL>
    if (~isOpen)
        mycluster = parcluster;
        matlabpool('local',mycluster.NumWorkers);
    end
end

fid = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          N OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'N')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the N Operator:   %dx%dx%d\n\n',Lr,Mr,Nr);
    

    [opToeplitzPlus,opToeplitzMinus] = assemblyC_Nop(r,r_offset,ko,dx);
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulantC_Nop(opToeplitzPlus,opToeplitzMinus); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = [opToeplitzPlus,opToeplitzMinus];           
    elseif strcmp(type ,'C')
        op_out = circulantC_Nop(opToeplitzPlus,opToeplitzMinus);    
    end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          K OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'K')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the K Operator:   %dx%dx%d\n\n',Lr,Mr,Nr);
    
    [opToeplitzPlus,opToeplitzMinus] = assemblyC_Kop(r,r_offset,ko,dx);
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulantC_Kop(opToeplitzPlus,opToeplitzMinus); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = [opToeplitzPlus,opToeplitzMinus];           
    elseif strcmp(type ,'C')
        op_out = circulant_Kop(opToeplitzPlus,opToeplitzMinus);    
    end
    
end

infocir = whos('op_out');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n     Dimensions:       %dx%dx%d',size(op_out,1),size(op_out,2), size(op_out,3));
fprintf(fid, '\n     Memory:           %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');

