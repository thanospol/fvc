function [J] = nearby_correct(m,I_non_correct,ST,EAc,EAo,VAc,VAo)
%%
J = zeros(6,6);
%
if isequal(m,[1,1,1]) %m == [1,1,1]; 
    
    C = [ ST  0.0 EAo EAo EAo EAo
          0.0 ST  EAo EAo EAo EAo
          EAo EAo ST  0.0 EAo EAo
          EAo EAo 0.0 ST  EAo EAo
          EAo EAo EAo EAo ST  0.0
          EAo EAo EAo EAo 0.0 ST  ];

    P = [ 0.0 1.0 0.0 0.0 0.0 0.0
          1.0 0.0 0.0 0.0 0.0 0.0
          0.0 0.0 0.0 1.0 0.0 0.0
          0.0 0.0 1.0 0.0 0.0 0.0
          0.0 0.0 0.0 0.0 0.0 1.0
          0.0 0.0 0.0 0.0 1.0 0.0 ];
      
elseif isequal(m,[2,1,1]) %m == [2,1,1]; 
    
    C = [ 0.0 ST  EAo EAo EAo EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 EAo EAc 0.0 VAo VAo
          0.0 EAo 0.0 EAc VAo VAo
          0.0 EAo VAo VAo EAc 0.0
          0.0 EAo VAo VAo 0.0 EAc ];

    P = [ 1.0 0.0 0.0 0.0 0.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 0.0 1.0 0.0 0.0
          1.0 0.0 1.0 0.0 0.0 0.0
          1.0 0.0 0.0 0.0 0.0 1.0
          1.0 0.0 0.0 0.0 1.0 0.0 ];
      
 elseif isequal(m,[2,1,2]) %m == [2,1,2]; 
    
    C = [ 0.0 EAc VAo VAo 0.0 EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo VAc 0.0 0.0 VAo
          0.0 VAo 0.0 VAc 0.0 VAo
          0.0 EAo VAo VAo 0.0 EAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];

    P = [ 1.0 0.0 0.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 0.0 1.0 1.0 0.0
          1.0 0.0 1.0 0.0 1.0 0.0
          1.0 0.0 0.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0 ];  
      
 elseif isequal(m,[1,1,2]) %m == [1,1,2]; 
    
    C = [ EAc 0.0 VAo VAo 0.0 EAo
          0.0 EAc VAo VAo 0.0 EAo
          VAo VAo EAc 0.0 0.0 EAo
          VAo VAo 0.0 EAc 0.0 EAo
          EAo EAo EAo EAo 0.0 ST
          0.0 0.0 0.0 0.0 0.0 0.0 ];

    P = [ 0.0 1.0 0.0 0.0 1.0 0.0
          1.0 0.0 0.0 0.0 1.0 0.0
          0.0 0.0 0.0 1.0 1.0 0.0
          0.0 0.0 1.0 0.0 1.0 0.0
          0.0 0.0 0.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0 ];  
      
 elseif isequal(m,[1,2,1]) %m == [1,2,1]; 
    
    C = [ EAc 0.0 0.0 EAo VAo VAo
          0.0 EAc 0.0 EAo VAo VAo
          EAo EAo 0.0 ST  EAo EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          VAo VAo 0.0 EAo EAc 0.0
          VAo VAo 0.0 EAo 0.0 EAc ];

    P = [ 0.0 1.0 1.0 0.0 0.0 0.0
          1.0 0.0 1.0 0.0 0.0 0.0
          0.0 0.0 1.0 0.0 0.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          0.0 0.0 1.0 0.0 0.0 1.0
          0.0 0.0 1.0 0.0 1.0 0.0 ]; 

 elseif isequal(m,[2,2,1]) %m == [2,2,1]; 
    
    C = [ 0.0 EAc 0.0 EAo VAo VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 EAo 0.0 EAc VAo VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAo VAc 0.0
          0.0 VAo 0.0 VAo 0.0 VAc ];

    P = [ 1.0 0.0 1.0 0.0 0.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 1.0 0.0 0.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 1.0 0.0 0.0 1.0
          1.0 0.0 1.0 0.0 1.0 0.0 ];  

 elseif isequal(m,[1,2,2]) %m == [1,2,2]; 
    
    C = [ VAc 0.0 0.0 VAo 0.0 VAo
          0.0 VAc 0.0 VAo 0.0 VAo
          VAo VAo 0.0 EAc 0.0 EAo
          0.0 0.0 0.0 0.0 0.0 0.0
          VAo VAo 0.0 EAo 0.0 EAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];

    P = [ 0.0 1.0 1.0 0.0 1.0 0.0
          1.0 0.0 1.0 0.0 1.0 0.0
          0.0 0.0 1.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          0.0 0.0 1.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0 ];  

 elseif isequal(m,[2,2,2]) %m == [2,2,2]; 
    
    C = [ 0.0 VAc 0.0 VAo 0.0 VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAc 0.0 VAo
          0.0 0.0 0.0 0.0 0.0 0.0
          0.0 VAo 0.0 VAo 0.0 VAc
          0.0 0.0 0.0 0.0 0.0 0.0 ];

    P = [ 1.0 0.0 1.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 1.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0
          1.0 0.0 1.0 0.0 1.0 0.0
          1.0 1.0 1.0 1.0 1.0 1.0 ]; 
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 6
    for jj = 1 : 6
        if P(ii,jj) == 0.0
            J(ii,jj) = C(ii,jj);
        elseif P(ii,jj) == 1.0
            J(ii,jj) = I_non_correct(ii,jj);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%