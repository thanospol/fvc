function absorb = j_inc(Etot,er,r,idx,freq,Nphoton)
%%  Supplies medium response to the pump for the N-photon process

%---local field intensity----%
[L,M,N,~] = size(r);
Length=L*M*N;

Eall=Etot(:);
Ex=Eall(1:Length);
Ey=Eall(Length+1:2*Length);
Ez=Eall(2*Length+1:3*Length);

Intensity=abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
%---absorbed function----%
jinc=-Intensity.^Nphoton.*imag(er(:)); %for N-photon process
d_jinc = [jinc;jinc;jinc];

jincidx = d_jinc(idx);
absorb    = spdiags( jincidx, 0, length(idx), length(idx) );

end