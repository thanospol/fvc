function epsilon = GST_permittivity(f,T)
% Dispersion relation for G2S2T5, as a function of frequency (Hz) and temperature (Kelvin).
% Data are from  Japanese Journal of Applied Physics 47, 5477 (2008), and Nature materials 7, 653 (2008).

x0=f/2.418e14; %unit to eV

Tsize=size(T);
epsilon=ones(Tsize);
T0=T-273.15;

Tamor=[26.85;163;];
Tcry=[191;348;]; % the temperature at which the  dispersion data is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load reamor1.csv;DATA=reamor1;
xn0=DATA(:,1);yn0=DATA(:,2);
xn=[xn0];yn=[yn0];
fe1=fit(xn,yn,'linearinterp');

load imamor.csv;DATA=imamor;
xk0=DATA(:,1);yk0=DATA(:,2);
xk=[0;xk0];yk=[0;yk0];
fe2=fit(xk,yk,'linearinterp');

er20=fe1(x0);
ep20=fe2(x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load n163.csv;DATA=n163;
xn0=DATA(:,1);yn0=DATA(:,2);
xn=[xn0];yn=[yn0];
fn1=fit(xn,yn,'linearinterp');

load k163.csv;DATA=k163;
xk0=DATA(:,1);yk0=DATA(:,2);
xk=[0;xk0];yk=[0;yk0];
fk1=fit(xk,yk,'linearinterp');

k0=fk1(x0);n0=fn1(x0);
ep=(n0+i*k0).^2;
er163=real(ep);ep163=imag(ep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load n191.csv;DATA=n191;
xn0=DATA(:,1);yn0=DATA(:,2);
xn=[xn0];yn=[yn0];
fn2=fit(xn,yn,'linearinterp');

load k191.csv;DATA=k191;
xk0=DATA(:,1);yk0=DATA(:,2);
xk=[0;xk0];yk=[0;yk0];
fk2=fit(xk,yk,'linearinterp');

k0=fk2(x0);n0=fn2(x0);
ep=(n0+i*k0).^2;
er191=real(ep);ep191=imag(ep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load n348.csv;DATA=n348;
xn0=DATA(:,1);yn0=DATA(:,2);
xn=[xn0];yn=[yn0];
fn3=fit(xn,yn,'linearinterp');

load k348.csv;DATA=k348;
xk0=DATA(:,1);yk0=DATA(:,2);
xk=[0;xk0];yk=[0;yk0];
fk3=fit(xk,yk,'linearinterp');

k0=fk3(x0);n0=fn3(x0);
ep=(n0+i*k0).^2;
er348=real(ep);ep348=imag(ep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load n440.csv;DATA=n440;
xn0=DATA(:,1);yn0=DATA(:,2);
xn=[xn0];yn=[yn0];
fn4=fit(xn,yn,'linearinterp');

load k440.csv;DATA=k440;
xk0=DATA(:,1);yk0=DATA(:,2);
xk=[0;xk0];yk=[0;yk0];
fk4=fit(xk,yk,'linearinterp');

k0=fk4(x0);n0=fn4(x0);
ep=(n0+i*k0).^2;
er440=real(ep);ep440=imag(ep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eramor=[er20;er163;];ercry=[er191;er348;];
epamor=[ep20;ep163;];epcry=[ep191;ep348;];

feramor=fit(Tamor,eramor,'linearinterp');
fepamor=fit(Tamor,epamor,'linearinterp');
fercry=fit(Tcry,ercry,'linearinterp');
fepcry=fit(Tcry,epcry,'linearinterp');

amphase=T0<165; % T=165 c is the phase transtion temperature from amorphour to cubic
tem=T0(amphase);
fcc=T0<350 & T0>=165;
temb=T0(fcc);
hex=T0>=350;% T=350 c is the phase transtion temperature from cubic to hexigonal
temh=T0(hex);

e1amor=feramor(tem);
e2amor=fepamor(tem);
epsilon(amphase)=e1amor-i*e2amor;

e1cry=fercry(temb);
e2cry=fepcry(temb);
epsilon(fcc)=e1cry-i*e2cry;
epsilon(hex)=er440-i*ep440;

