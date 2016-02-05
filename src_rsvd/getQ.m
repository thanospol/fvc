function [Q,ii] = getQ(Fdirect,m,n,Lmax,tol,blocksize)
%% Compute the approximate QR decomposition

t2 = tic;

fid = 1;
% if on-the-fly control of SV drop
if (Lmax < 1)
    
    % allocate memory for the blocksize
    AGn = ones(m,blocksize) + 1j*ones(m,blocksize);
    Q = [];
    
    flagdone = 0;
    jj = 0;
    ii = 0;    
    while (flagdone == 0)
        
        An_v = randn(n,1) + 1j*randn(n,1);
        An_v= An_v/norm(An_v); % normalize excitation
                        
        % compute direct operator
        jj = jj + 1;
        ii = ii + 1;
        [AGn(:,jj)] = Fdirect(An_v);
        
        if ((ii/blocksize) == floor(ii/blocksize))
            
            % optional on-the-fly control of SV drop
            
            if size(Q,2) % there are already vectors in Q
                AGn = AGn - Q*(Q'*AGn); % orthogonalize newblock w.r.t. Q
                [AGn, Sr, ~] = svd(AGn, 'econ'); % apply SVD to newblock
                Q = [Q, AGn]; % add vectors to Q
                
            else
                [Q, Sr, ~] = svd(AGn, 'econ'); % apply SVD to newblock
                Smax = Sr(1,1);
            end
            
            % check if newblock is rank deficient for the defined tolerance
            flagdone = 0;
            for kk = 1:length(Sr)
                if Sr(kk,kk) < 0.5*tol*Smax % use looser tolerance
                    flagdone = 1; % if so, we are done
                    break
                end
            end
            
            jj = 0;
            fprintf(fid, '\n %4.0d random excitations done (Max SV = %g, current SV drop ratio %1.12f). Elapsed time %g', ii, Smax, Sr(kk,kk)/Smax,toc(t2));
            t2 = tic;
            
        end
       
    end
    clear Einc;

    
else % no control... predefined set of excitations Lmax
    
    % allocate memory
    An = randn(n,Lmax) + 1j*randn(n,Lmax);
    AGn = ones(m,Lmax) + 1j*ones(m,Lmax);
    
    for ii = 1:Lmax
        
        An_v = An(:,ii)/norm(An(:,ii)); % normalize excitation
        An(:,ii) = An_v;
        % compute direct operator
        [AGn(:,ii)] = Fdirect(An_v);
        
        if ((ii/blocksize) == floor(ii/blocksize))
            fprintf(fid, '\n %4.0d random excitations done. Elapsed time %g', ii, toc(t2));
            t2 = tic;
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Compute the SVD of the basis (equivalent to QR)
    [Q, ~, ~] = svd(AGn,'econ');
%     [Q, ~] = qr(AGn,0);
   
end
