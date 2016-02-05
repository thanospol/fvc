
% Nop
mex -output mexCUBATURE_Nop -v -O Cubature_Nop_mex.cpp Cubature_Nop.cpp
% Kop
mex -output mexCUBATURE_Kop -v -O Cubature_Kop_mex.cpp Cubature_Kop.cpp 
% DIRECTFN
mex -output mexDIRECTFN_WS_ST_const -largeArrayDims -v -O DIRECTFN_WS_ST_const_mex.cpp DIRECTFN_WS_ST_const.cpp
mex -output mexDIRECTFN_WS_EA_const -largeArrayDims -v -O DIRECTFN_WS_EA_const_mex.cpp DIRECTFN_WS_EA_const.cpp 
mex -output mexDIRECTFN_WS_VA_const -largeArrayDims -v -O DIRECTFN_WS_VA_const_mex.cpp DIRECTFN_WS_VA_const.cpp 