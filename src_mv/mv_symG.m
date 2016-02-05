function [JOut] = mv_symG(JIn, fN, idx, ce)
%% Function that applies the operator symG: (G + ctranspose(G))/2

% the non-transpose component
[J_nt] = mv_AN(JIn, fN, 0, -1, 0, ce, idx, 'notransp');

% the transpose component
[J_t] = mv_AN(JIn, fN, 0, -1, 0, ce, idx, 'transp');

% Output
JOut = (  J_nt +  J_t)/2;

end
