%% instal fvc


%% mex functions
%
% First, make sure you have assigned a c compiler for matlab: mex -setup
%

ST = 'mexDIRECTFN_WS_ST_const';
EA = 'mexDIRECTFN_WS_EA_const';
VA = 'mexDIRECTFN_WS_VA_const';
Nop = 'mexCUBATURE_Nop';
Kop = 'mexCUBATURE_Kop';

if (exist(ST,'file')==3 && exist(EA,'file')==3 && exist(VA,'file')==3 && exist(Nop,'file')==3 && exist(Kop,'file')==3)
    
    display('mex functions already exist')
    
else % src_operators\src_mexcubatures
    
    if ispc
        
        currentfolder = pwd;
        
        cd .\src_operators\src_mexcubatures
        
        make_fvc
        
        cd(currentfolder)
        
    else
        
        currentfolder = pwd;
        
        cd .//src_operators/src_mexcubatures
        
        make_fvc
        
        cd(currentfolder)
    end
    
end


%% save path

p = genpath(pwd);
addpath(p);
savepath;

