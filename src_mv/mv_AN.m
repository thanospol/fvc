function [JOut] = mv_AN(JIn0, fN, Mr, Mc, Gram, ce, idx, transp_flag)
%% Function that applies the operator AN: (1/ce) * ( Mr*Gram - Mc*N )

% fft dimensions
LfN = size(fN,1);
MfN = size(fN,2);
NfN = size(fN,3);

L = LfN/2;
M = MfN/2;
N = NfN/2;

% allocate space
JIn  = zeros(L, M, N, 3);
JOut = zeros(L, M, N, 3);

if(nargin == 6 || (isempty(idx) && isempty(transp_flag)) )
   idx = 1:(3*L*M*N);
   transp_flag = 'notransp';
end

if(nargin == 7 || isempty(transp_flag) )
   transp_flag = 'notransp';
end

% translate from local to global coordinates
JIn(idx) = JIn0(:);

if strcmp(transp_flag,'notransp') % y = Mchi*N*x
    
    % apply fft and mv-op
    fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
    
    % apply ifft and form the final product
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = Gram .* Mr .* JIn(:,:,:,1) - Mc .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = Gram .* Mr .* JIn(:,:,:,2) - Mc .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = Gram .* Mr .* JIn(:,:,:,3) - Mc .* Jout3(1:L,1:M,1:N);
    
    % Output
    JOut = (1/ce) * JOut;
    
elseif strcmp(transp_flag,'conjugate') % y = N * conj(Mchi)*x 
    
    % apply fft and mv-op
    fJ = fftn( Mc .* conj(JIn(:,:,:,1)),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;
    
    fJ = fftn( Mc .* conj(JIn(:,:,:,2)),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    fJ = fftn( Mc .* conj(JIn(:,:,:,3)),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
       
    % apply ifft and form the final product
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = Gram .* conj(Mr) .* JIn(:,:,:,1) - conj( Jout1(1:L,1:M,1:N) );
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = Gram .* conj(Mr) .* JIn(:,:,:,2) - conj( Jout2(1:L,1:M,1:N) );
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = Gram .* conj(Mr) .* JIn(:,:,:,3) - conj( Jout3(1:L,1:M,1:N) );
    
    % Output
    JOut = conj(1/ce) * JOut;

elseif strcmp(transp_flag,'conjugate2') % y = conj(Mchi*N)*x 
    
    % apply fft and mv-op
    fJ = fftn(conj(JIn(:,:,:,1)),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;
    
    fJ = fftn(conj(JIn(:,:,:,2)),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    fJ = fftn(conj(JIn(:,:,:,3)),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
       
    % apply ifft and form the final product
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = Gram .* conj(Mr) .* JIn(:,:,:,1) - conj(Mc .* Jout1(1:L,1:M,1:N) );
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = Gram .* conj(Mr) .* JIn(:,:,:,2) - conj(Mc .* Jout2(1:L,1:M,1:N) );
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = Gram .* conj(Mr) .* JIn(:,:,:,3) - conj(Mc .* Jout3(1:L,1:M,1:N) );
    
    % Output
    JOut = conj(1/ce) * JOut;
       
elseif strcmp(transp_flag,'transp_unconjugated') % y = A.'*x = flip(A*flip(x))
    
    % -------------------------------------------------------------------------
    % first compute transpose(Mc*N)*JIn
    % -------------------------------------------------------------------------
    
    JIn2(:,:,:,1) = Mc.* JIn(:,:,:,1);
    JIn2(:,:,:,2) = Mc.* JIn(:,:,:,2);
    JIn2(:,:,:,3) = Mc.* JIn(:,:,:,3);
    
    % flip J
    JIn2 = flipJ(JIn2,L,M,N);
    
    % apply fft and mv-op
    fJ = fftn(JIn2(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn2(:,:,:,2),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn2(:,:,:,3),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = Jout3(1:L,1:M,1:N);
    
    % flip Jout
    JOut = flipJ(JOut,L,M,N);
    
    % -------------------------------------------------------------------------
    % transpose(Mc*Gram)* JIn - transpose(Mc*N)*JIn
    % -------------------------------------------------------------------------
    JOut(:,:,:,1) = Gram .* Mr .* JIn(:,:,:,1) - JOut(:,:,:,1);
    JOut(:,:,:,2) = Gram .* Mr .* JIn(:,:,:,2) - JOut(:,:,:,2);
    JOut(:,:,:,3) = Gram .* Mr .* JIn(:,:,:,3) - JOut(:,:,:,3);
    
    % Output
    JOut = (1/ce) * JOut;
    
elseif strcmp(transp_flag,'transp') % y = A'*x = conj(A.'*conj(x)) = conj(flip(A*flip(conj(x))))

    
    % -------------------------------------------------------------------------
    % first compute (Mc*N)' * JIn = N' * (Mc'* JIn)
    % -------------------------------------------------------------------------
    
    JIn2(:,:,:,1) = conj(Mc).* JIn(:,:,:,1);
    JIn2(:,:,:,2) = conj(Mc).* JIn(:,:,:,2);
    JIn2(:,:,:,3) = conj(Mc).* JIn(:,:,:,3);
    
    % flip conj(J)
    JIn2 = flipJ(conj(JIn2),L,M,N);
    
    % apply fft and mv-op
    fJ = fftn(JIn2(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = fN(:,:,:,1) .* fJ;
    Jout2 = fN(:,:,:,2) .* fJ;
    Jout3 = fN(:,:,:,3) .* fJ;
    
    fJ = fftn(JIn2(:,:,:,2),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,2) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,4) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,5) .* fJ;
    
    fJ = fftn(JIn2(:,:,:,3),[LfN, MfN, NfN]);
    Jout1 = Jout1 + fN(:,:,:,3) .* fJ;
    Jout2 = Jout2 + fN(:,:,:,5) .* fJ;
    Jout3 = Jout3 + fN(:,:,:,6) .* fJ;
    
    % apply ifft
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = conj(Jout1(1:L,1:M,1:N)) ;
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = conj(Jout2(1:L,1:M,1:N)) ;
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = conj(Jout3(1:L,1:M,1:N)) ;
    
    % flip Jout
    JOut = flipJ(JOut,L,M,N);
    
    % -------------------------------------------------------------------------
    % (Mc*Gram)' * JIn - (Mc*N)' * JIn
    % -------------------------------------------------------------------------
    JOut(:,:,:,1) = Gram .* conj(Mr) .* JIn(:,:,:,1) - JOut(:,:,:,1);
    JOut(:,:,:,2) = Gram .* conj(Mr) .* JIn(:,:,:,2) - JOut(:,:,:,2);
    JOut(:,:,:,3) = Gram .* conj(Mr) .* JIn(:,:,:,3) - JOut(:,:,:,3);
    
    % Output
    JOut = conj(1/ce) * JOut;
    
end

% Output
JOut = reshape(JOut,3*L*M*N,1);
JOut = JOut(idx);

end

%% local function
function Jf = flipJ(J,L,M,N)
    
    Js = J(:,:,:,1);
    Js = Js(end:-1:1); % Js = flip(Js(:)); for Matlab > 2014
    Jf(:,:,:,1) = reshape(Js,L,M,N);
    
    Js = J(:,:,:,2);
    Js = Js(end:-1:1); % Js = flip(Js(:)); for Matlab > 2014
    Jf(:,:,:,2) = reshape(Js,L,M,N);
    
    Js = J(:,:,:,3);
    Js = Js(end:-1:1); % Js = flip(Js(:)); for Matlab > 2014
    Jf(:,:,:,3) = reshape(Js,L,M,N);
    

end
