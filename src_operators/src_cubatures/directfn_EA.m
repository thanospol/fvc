function [I_EAc,I_EAo] = directfn_EA(Np_1D,ko,dx,dy,dz)
%%

% GEOMETRY
r1 = [0, 0,0]';
r2 = [dx,0,0]';
r3 = [dx,dy,0]';
r4 = [0 ,dy,0]';
r5 = [2*dx,0,0]';
r6 = [2*dx ,dy,0]';
r5o = [dx,0,dz]';
r6o = [dx ,dy,dz]';
%
% N_theta = Np_1D;
% N_psi = Np_1D;
% %
% [ w_theta , z_theta ] = Gauss_1D ( N_theta );
% [ w_psi , z_psi ] = Gauss_1D ( N_psi );
% %
% N_theta_p   = Np_1D;
% N_theta_q   = Np_1D;
% N_psi       = Np_1D;
% %
% [ w_theta_p,z_theta_p ] = Gauss_1D ( N_theta_p );
% [ w_theta_q,z_theta_q ] = Gauss_1D ( N_theta_q );
% [ w_psi,z_psi ]         = Gauss_1D_Nop ( N_psi );
% Evaluate EAc

[I_T1_T1] = mexDIRECTFN_WS_EA_const(r2,r3,r1,r5, ko, Np_1D);
[I_T1_T2] = mexDIRECTFN_WS_VA_const(r3,r1,r2,r5,r6, ko, Np_1D);
[I_T2_T1] = mexDIRECTFN_WS_VA_const(r3,r4,r1,r2,r5, ko, Np_1D);
[I_T2_T2] = mexDIRECTFN_WS_VA_const(r3,r4,r1,r5,r6, ko, Np_1D);

%
I_EAc = I_T1_T1 + I_T1_T2 + I_T2_T1 + I_T2_T2;

% Evaluate EAo
[I_T1_T1] = mexDIRECTFN_WS_EA_const(r2,r3,r1,r5o, ko, Np_1D);
[I_T1_T2] = mexDIRECTFN_WS_VA_const(r3,r1,r2,r5o,r6o, ko, Np_1D);
[I_T2_T1] = mexDIRECTFN_WS_VA_const(r3,r4,r1,r2,r5o, ko, Np_1D);
[I_T2_T2] = mexDIRECTFN_WS_VA_const(r3,r4,r1,r5o,r6o, ko, Np_1D);

%
I_EAo = I_T1_T1 + I_T1_T2 + I_T2_T1 + I_T2_T2;