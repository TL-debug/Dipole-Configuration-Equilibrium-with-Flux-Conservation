function [GBR,GBZ,DGBR,DGBZ] = GradGrep(Coils,sep,dr,dz)
% GRADGREP 对固定点和线圈的格林函数求一阶和二阶导数
% 输入: 线圈位置，固定点位置，横向和垂直步长
% 输出: 格林函数的一阶和二阶导数
% Note: 线圈的位置是固定的，选择固定点的周围点进行求导
%       我们采用先计算计算域的格林函数的导数，在插值给出固定点的各个量

% GBR = 0*R;GBZ = GBR;DGBR = GBR;DGBZ = GBR;
% NR = size(R,1);NZ = size(R,2);
% GrP = Grep(Coils,sep);
% 
% for i = 1:NR
%     for j = 1:NZ
% 
%         if ( j>=3 && j<=NZ-2 && i>=3 && i<=NR-2  )
% 
%             GBR(i,j)  = -(-1/12*GrP(i,j+2) + 8/12*GrP(i,j+1) - ...
%                           8/12*GrP(i,j-1) + 1/12*GrP(i,j-2))/(R(i,j)*dz);
%             GBZ(i,j)  =  (-1/12*GrP(i+2,j) + 8/12*GrP(i+1,j) - ...
%                           8/12*GrP(i-1,j) + 1/12*GrP(i-2,j))/(R(i,j)*dr);
%             DGBR(i,j) = -(-1/12*GrP(i,j+2) + 16/12*GrP(i,j+1) - 30/12*GrP(i,j) +...
%                           16/12*GrP(i,j-1) - 1/12*GrP(i,j-2))/(R(i,j)*dz^2);
%             DGBZ(i,j) =  (-1/12*GrP(i+2,j) + 16/12*GrP(i+1,j) - 30/12*GrP(i,j) +...
%                           16/12*GrP(i-1,j) - 1/12*GrP(i-2,j))/(R(i,j)*dr^2);
%         else 
%             GBR(i,j)  = -(GrP(i,j+1) - GrP(i,j-1))/(2*R(i,j)*dz);
%             GBZ(i,j)  =  (GrP(i+1,j) - GrP(i-1,j))/(2*R(i,j)*dr);
%             DGBR(i,j) = -(GrP(i,j+1) - 2*GrP(i,j) + GrP(i,j-1))/(R(i,j)*dz^2);
%             DGBZ(i,j) =  (GrP(i+1,j) - 2*GrP(i,j) + GrP(i-1,j))/(R(i,j)*dr^2);
%         end
%     end
% end
% 
% tGBR  = interp2(R',Z',GBR',)
% tGBZ  =
% tDGBR =
% tDGBZ = 




% - 一阶导数：Grep返回行表示固定点，列对应线圈
Gr = Grep(Coils,sep);
tsep = sep;tsep(1,:) = tsep(1,:) + 2*dr;Gppr = Grep(Coils,tsep);
tsep = sep;tsep(1,:) = tsep(1,:) + dr;  Gpr  = Grep(Coils,tsep);
tsep = sep;tsep(1,:) = tsep(1,:) - dr;  Gmr  = Grep(Coils,tsep);
tsep = sep;tsep(1,:) = tsep(1,:) - 2*dr;Gmmr = Grep(Coils,tsep);

tsep = sep;tsep(2,:) = tsep(2,:) + 2*dz;Gppz = Grep(Coils,tsep);
tsep = sep;tsep(2,:) = tsep(2,:) + dz;  Gpz  = Grep(Coils,tsep);
tsep = sep;tsep(2,:) = tsep(2,:) - dz;  Gmz  = Grep(Coils,tsep);
tsep = sep;tsep(2,:) = tsep(2,:) - 2*dz;Gmmz = Grep(Coils,tsep);

MRR = repmat(sep(1,:)',1,size(Coils,1));

GBR = -(-1/12*Gppz + 8/12*Gpz - 8/12*Gmz + 1/12*Gmmz)./(MRR*dz);
GBZ =  (-1/12*Gppr + 8/12*Gpr - 8/12*Gmr + 1/12*Gmmr)./(MRR*dr);

% - 二阶导数
DGBR =  (1/12*Gppz - 16/12*Gpz + 30/12*Gr - 16/12*Gmz + 1/12*Gmmz)./(MRR*dz^2);
DGBZ = -(1/12*Gppr - 16/12*Gpr + 30/12*Gr - 16/12*Gmr + 1/12*Gmmr)./(MRR*dr^2)-...
        (-1/12*Gppr + 8/12*Gpr - 8/12*Gmr + 1/12*Gmmr)./(MRR*dr);% 这里与原文中的导数存在差异



end

