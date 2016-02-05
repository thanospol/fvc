function plotGeometry(r,EMT)
%% plot the geometry under study
idx = find(EMT.Er-1);

x = r(:,:,:,1);
y = r(:,:,:,2);
z = r(:,:,:,3);

% figure
figure()
plot3(x(idx),y(idx),z(idx),'bo');
daspect([1 1 1])
grid on

end