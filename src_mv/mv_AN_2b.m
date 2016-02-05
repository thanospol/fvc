function [y] = mv_AN_2b(x, fN_11, fN_12, EMT_1, EMT_2, M1, M2, Gram, transp_flag)
%% Function that applies the operator AN for 2 bodies
%
% INPUT
% JIn:  Current 
% fN:   FFT-Circulant of N operator
% Mer:  Relative lectric permittivity  
% Mer:  Relative magnetic permeability 
% Mce:  Electric susceptibility 
% Mcm:  Magnetic susceptibility 
% Gram: Gram matrix 
% ce:   constant jwe
% ce:   constant jwm 
%
% transp_flag: flag for choosing the m-v product
%         'transp':   A' * b
%         'notransp': A * b
%
% idx: index with local coordinates (non-air voxels)


% domain dimensions
% [L1, M1, N1] = size(EMT_1.Er);
% [L2, M2, N2] = size(EMT_2.Er);
% 
% sizeJ = 3*L*M*N;
% 
% Je = JIn(1:sizeJ,1);
% Jm = JIn((sizeJ+1):end,1);


% N22 = N11
% N21 = transpose(N12)

ce = 1;

sz_x1 = length(M1.idx);
sz_x2 = length(M2.idx);

x1 = x(1:sz_x1,1);
x2(1:sz_x2,1) = x(sz_x1+1:end,1);

if strcmp(transp_flag,'transp')   

    [y1] = mv_AN(x1, fN_11, EMT_1.Er, EMT_1.Er-1, Gram, ce, M1.idx, 'transp') +  mv_AN(x2, fN_12, EMT_2.Er, EMT_2.Er-1, 0   , ce, M2.idx, 'conjugate');
    [y2] = mv_AN(x1, fN_12, EMT_1.Er, EMT_1.Er-1, 0   , ce, M1.idx, 'transp') +  mv_AN(x2, fN_11, EMT_2.Er, EMT_2.Er-1, Gram, ce, M2.idx, 'transp');

elseif strcmp(transp_flag,'notransp')
 
    [y1] = mv_AN(x1, fN_11, EMT_1.Er, EMT_1.Er-1, Gram, ce, M1.idx, 'notransp') +  mv_AN(x2, fN_12, EMT_1.Er, EMT_1.Er-1, 0   , ce, M2.idx, 'notransp');
    [y2] = mv_AN(x1, fN_12, EMT_2.Er, EMT_2.Er-1, 0   , ce, M1.idx, 'transp_unconjugated')   +  mv_AN(x2, fN_11, EMT_2.Er, EMT_2.Er-1, Gram, ce, M2.idx, 'notransp');
    
end

y = [y1;y2];

end


