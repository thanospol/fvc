function [G_mn1,G_mn2] = assemblyC_Nop(r,r_offset,ko,dx)
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Grid Dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, M, N, ~] = size(r);
dy =dx; dz = dx;

% Get the 6 faces of  the cube
[R_faces] = cube_faces(dx,dy,dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np_1D_far    = 4;

% Allocate memory for main matrices
J_mn1 = zeros(L,M,N,6,6);
G_mn1 = zeros(L,M,N,6);

J_mn2 = zeros(L,M,N,6,6);
G_mn2 = zeros(L,M,N,6);
% Reference cell
r_n = [0.0 , 0.0 , 0.0]';
% Reference distance vector
d = [dx,dy,dz];

tic_Assembly = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_far;

[Np,wp,up,vp] = gauss_2D(N_P);

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Far Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = R_faces;
kko = ko;
ddx = dx;
%
parfor mx = 1:L
    for my = 1:M
        for mz = 1:N
            m = [mx,my,mz];
            r_m_plus  = ( (m-1) .* d + r_offset )';
            r_m_minus = ( (m-1) .* d - r_offset )';
                        

            J_mn1(mx,my,mz,:,:) = mexCUBATURE_Nop(Np,wp,up,vp,r_m_plus,r_n,RR,kko,ddx);
            J_mn2(mx,my,mz,:,:) = mexCUBATURE_Nop(Np,wp,up,vp,r_m_minus,r_n,RR,kko,ddx);
                     
        end
    end
end

Time_far = toc;
fprintf('Time_far         = %d \n',int64(Time_far));
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Evaluate G_mn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
tic

[P_vector] = vector_projection();

for kk = 1:6
    for jj = 1:6
        for ii = 1:6
            G_mn1(:,:,:,kk) = G_mn1(:,:,:,kk) + J_mn1(:,:,:,ii,jj)*P_vector(ii,jj,kk);
            G_mn2(:,:,:,kk) = G_mn2(:,:,:,kk) + J_mn2(:,:,:,ii,jj)*P_vector(ii,jj,kk);
        end
    end
end

Time_projection = toc;
fprintf('Time_projection  = %d \n',int64(Time_projection));

G_mn1 = ( 1.0 / (4*pi) ) * G_mn1;
G_mn2 = ( 1.0 / (4*pi) ) * G_mn2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_Assembly = toc(tic_Assembly);
fprintf('------------------------  \n')
fprintf('Time_Assembly    = %dm%ds  \n' ,floor(Time_Assembly/60),int64(mod(Time_Assembly,60)))
