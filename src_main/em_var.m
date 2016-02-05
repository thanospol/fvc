function EM =  em_var(f)
%% some basic electromagnetic constansts and variables

EM.mo = 4*pi*1e-7;
EM.co = 299792458;
EM.eo = 1/EM.co^2/EM.mo;
%
EM.omega = 2 * pi * f;
EM.lambda  = EM.co/f;
EM.ko = 2*pi/EM.lambda;
EM.omega_mu = EM.omega * EM.mo;
EM.eta =  3.767303134617706e+002; 
%
EM.ce = 1i*EM.omega*EM.eo;
EM.cm = 1i*EM.omega*EM.mo;

end