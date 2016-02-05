function plot_LDOS(r,er,Radius,LDOS)
%% plot the LDOS in X-Z cross-section

d_chiEx  = er(:) - 1.0;
idx0 = find(d_chiEx~=0);    
L=length(idx0);

X0=r(:,:,:,1);Y0=r(:,:,:,2);Z0=r(:,:,:,3);
X0=X0(:);Y0=Y0(:);Z0=Z0(:);
LDOS0=zeros(size(X0));

X=X0(idx0);Y=Y0(idx0);Z=Z0(idx0);
LDOS0(idx0)=LDOS(1:L)+LDOS(L+1:2*L)+LDOS(2*L+1:3*L); %adding the X,Y,Z components
F=scatteredInterpolant(X0,Y0,Z0,LDOS0,'linear');

Da=20;DD=Radius/Da;
XX=-Radius:DD:Radius;
[XX,ZZ]=meshgrid(XX,XX);YY=zeros(size(XX));

DF0=F(XX,YY,ZZ);DFM=max(DF0(:));
DF=DF0/DFM;%normalization

figure;
contourf(XX/1e-6,ZZ/1e-6,DF,'LineStyle','none');
colormap hot;
xlabel('X(\mu m)');ylabel('Z(\mu m)');
set(gcf, 'renderer', 'zbuffer');
colorbar
title('LDOS in X-Z cross section')
set(gca,'FontSize',20)
