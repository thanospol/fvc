function [r_1,r_2,EMT_1,EMT_2] = getGeometry_2obj(points_min,shape,Radius,r_offset)


% Object 1
Obj1Properties.Epsilon      = 12-1i*1;
Obj1Properties.Mu           = 1;
Obj1Properties.Temp_profile = 1;
Obj1Properties.shape        = shape;
Obj1Properties.Radius       = Radius;

% Object 2
Obj2Properties.Epsilon      = 12-1i*1;
Obj2Properties.Mu           = 1;
Obj2Properties.Temp_profile = 1;
Obj2Properties.shape        = shape;
Obj2Properties.Radius       = Radius;

% offset between the 2 objects
% r_offset = [2*micron+micron , 0.0, 0.0];

%
[r_1,EMT_1] = getGeometry_1obj(points_min,Obj1Properties);
[r_2,EMT_2] = getGeometry_1obj(points_min,Obj2Properties);
%
Roffset = ones(size(r_2));
Roffset(:,:,:,1) = Roffset(:,:,:,1)*r_offset(1);
Roffset(:,:,:,2) = Roffset(:,:,:,2)*r_offset(2);
Roffset(:,:,:,3) = Roffset(:,:,:,3)*r_offset(3);

r_2 = r_2 + Roffset;