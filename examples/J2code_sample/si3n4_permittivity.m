function epsilon = si3n4_permittivity(f)
% Dispersion relation for silicon nitride, as a function of frequency (Hz)
% Data are from  Applied optics 51, 6789 (2012).

load si3n4.txt;DATA=si3n4;
xn0=DATA(:,1);yn0=DATA(:,2);zn0=DATA(:,3);
fn0=3e8./(1e-6*xn0);
fn=[fn0;0];yn=[yn0;1];zn=[zn0;0];
nn=fit(fn,yn,'linearinterp');
kn=fit(fn,zn,'linearinterp');

n1=nn(f);
k1=kn(f);

epsilon=(n1-i*k1).^2;
