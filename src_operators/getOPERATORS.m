function [op_out] = getOPERATORS(r,f,op,type)
%% function that retrieves the operators N and K for JM-VIE
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
%
%% OUTPUT
% op_out    discrete operator
%
% -------------------------------------------------------------------------
%
%   A. G. Polimeridis -- thanos_p@mit.edu
%   J. Fernandez Villena -- jvillena@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 || isempty(type))
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
check_matlabpool;

fid = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          N OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'N')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the N Operator:   %dx%dx%d\n\n',Lr,Mr,Nr);
    
    [opToeplitz] = assembly_Nop(r,ko,dx);
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulant_Nop(opToeplitz); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = opToeplitz;           
    elseif strcmp(type ,'C')
        op_out = circulant_Nop(opToeplitz);    
    end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          K OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'K')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the K Operator:   %dx%dx%d\n\n',Lr,Mr,Nr);
    
    [opToeplitz] = assembly_Kop(r,ko,dx);
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulant_Kop(opToeplitz); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = opToeplitz;           
    elseif strcmp(type ,'C')
        op_out = circulant_Kop(opToeplitz);    
    end
    
end

infocir = whos('op_out');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n     Dimensions:       %dx%dx%d',size(op_out,1),size(op_out,2), size(op_out,3));
fprintf(fid, '\n     Memory:           %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');

