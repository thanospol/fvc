function [E] = e_field_Nop_comp(Jidx, fG,  dV, omega, eo, Einc, idx)
%% Compute the total electric fields due to electric currents
% E = Einc + Escat
%   = Einc + (1/(j*omega*eo)) N*J

[L, M, N, ~] = size(Einc);
J = zeros( L, M, N, 3);

J(idx) = Jidx(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = zeros(L, M, N, 3);

LfG = size(fG,1);
MfG = size(fG,2);
NfG = size(fG,3);


% Compute FFT(J) and mv-op
fJ = fftn(J(:,:,:,1), [LfG, MfG, NfG]);
Jout1 = fG(:,:,:,1) .* fJ;
Jout2 = fG(:,:,:,2) .* fJ;
Jout3 = fG(:,:,:,3) .* fJ;

fJ = fftn(J(:,:,:,2), [LfG, MfG, NfG]);
Jout1 = Jout1 + fG(:,:,:,2) .* fJ;
Jout2 = Jout2 + fG(:,:,:,4) .* fJ;
Jout3 = Jout3 + fG(:,:,:,5) .* fJ;

fJ = fftn(J(:,:,:,3), [LfG, MfG, NfG]);
Jout1 = Jout1 + fG(:,:,:,3) .* fJ;
Jout2 = Jout2 + fG(:,:,:,5) .* fJ;
Jout3 = Jout3 + fG(:,:,:,6) .* fJ;

% apply ifft
fJ = ifftn(Jout1);
E(:,:,:,1) = fJ(1:L, 1:M, 1:N);
fJ = ifftn(Jout2);
E(:,:,:,2) = fJ(1:L, 1:M, 1:N);
fJ = ifftn(Jout3);
E(:,:,:,3) = fJ(1:L, 1:M, 1:N);

% assembly field
E = E - dV * J;
E = (1/dV) * E;
E = (1/dV) * Einc + E./(1i*omega*eo);
end