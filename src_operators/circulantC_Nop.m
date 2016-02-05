function [C] = circulantC_Nop(Gp,Gm)
%%
[L, M, N, ~] = size(Gp);
%
C = zeros(2*L,2*M,2*N ,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cube 'L'
[Gp_coeff_L] = gperiodic_coeff_Nop('L');
% Cube 'M'
[Gp_coeff_M] = gperiodic_coeff_Nop('M');
% Cube 'N'
[Gp_coeff_N] = gperiodic_coeff_Nop('N');
% Cube 'LM'
[Gp_coeff_LM] = gperiodic_coeff_Nop('LM');
% Cube 'LN'
[Gp_coeff_LN] = gperiodic_coeff_Nop('LN');
% Cube 'MN'
[Gp_coeff_MN] = gperiodic_coeff_Nop('MN');
% Cube 'LMN' equal to 1
% [Gp_coeff_LMN] = gperiodic_coeff('LMN');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C(1:L,1:M,1:N,:) = Gp; 
%
for ii = 1:6
    % Cube 'L'
    C(L+2:2*L,1:M,1:N,ii)         = Gm(L:-1:2,1:M,1:N,ii) * Gp_coeff_L(ii,1);
    % Cube 'M'
    C(1:L,M+2:2*M,1:N,ii)         = Gp(1:L,M:-1:2,1:N,ii) * Gp_coeff_M(ii,1);
    % Cube 'N'
    C(1:L,1:M,N+2:2*N,ii)         = Gp(1:L,1:M,N:-1:2,ii) * Gp_coeff_N(ii,1);
    % Cube 'LM'
    C(L+2:2*L,M+2:2*M,1:N,ii)     = Gm(L:-1:2,M:-1:2,1:N,ii) * Gp_coeff_LM(ii,1);
    % Cube 'LN'
    C(L+2:2*L,1:M,N+2:2*N,ii)     = Gm(L:-1:2,1:M,N:-1:2,ii) * Gp_coeff_LN(ii,1);
    % Cube 'MN'
    C(1:L,M+2:2*M,N+2:2*N,ii)     = Gp(1:L,M:-1:2,N:-1:2,ii) * Gp_coeff_MN(ii,1);
end
% Cube 'LMN'
C(L+2:2*L,M+2:2*M,N+2:2*N,:) = Gm(L:-1:2,M:-1:2,N:-1:2,:) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%