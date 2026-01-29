function [FBR,FBZ,DFBR,DFBZ] = GradFB(sep,tmpsi,R,Z,dr,dz)
% GRADFB 固定点的磁场和磁场的导数
% 输入 固定点，通量，横向和纵向的网格，横向和纵向步长（sep只用到前两行）
% 输出 固定点的磁场横向和纵分量，一阶导数的横向和纵向分量
tmpBR = zeros(size(R));tmpBZ = tmpBR;
tmpdBR = tmpBR;tmpdBZ = tmpBZ;
NR = size(R,1);
NZ = size(R,2);
% 磁场分量和一阶导数
for i = 2:NR-1
    for j = 2:NZ-1
        if ( j>=3 && j<=NZ-2 && i>=3 && i<=NR-2 ) %五点差分
            % - 一阶导数
            tmpBR(i,j) = -(-1/12*tmpsi(i,j+2)+2/3*tmpsi(i,j+1)-2/3*tmpsi(i,j-1)+...
                          1/12*tmpsi(i,j-2))/(dz*R(i,j));
            tmpBZ(i,j) =  (-1/12*tmpsi(i+2,j)+2/3*tmpsi(i+1,j)-2/3*tmpsi(i-1,j)+...
                          1/12*tmpsi(i-2,j))/(dr*R(i,j));
            % - 二阶导数
            tmpdBR(i,j) =  (1/12*tmpsi(i,j+2)-16/12*tmpsi(i,j+1)+30/12*tmpsi(i,j)-...
                            16/12*tmpsi(i,j-1)+1/12*tmpsi(i,j-2))/(dz^2*R(i,j));
            tmpdBZ(i,j) = -(1/12*tmpsi(i+2,j)-16/12*tmpsi(i+1,j)+30/12*tmpsi(i,j)-...
                            16/12*tmpsi(i-1,j)+1/12*tmpsi(i-2,j))/(dr^2*R(i,j))-...
                           (-1/12*tmpsi(i+2,j)+2/3*tmpsi(i+1,j)-2/3*tmpsi(i-1,j)+...
                            1/12*tmpsi(i-2,j))/(dr*R(i,j)^2);
        elseif ( j>=2 && j<=NZ-1 && i>=2 && i<=NR-1 )
            tmpBR(i,j) = -(tmpsi(i,j+1)-tmpsi(i,j-1))/(2*dz*R(i,j));
            tmpBZ(i,j) =  (tmpsi(i+1,j)-tmpsi(i-1,j))/(2*dr*R(i,j));

            tmpdBR(i,j) = -(tmpsi(i,j+1)-2*tmpsi(i,j)+tmpsi(i,j-1))/(dz^2*R(i,j));
            tmpdBZ(i,j) =  (tmpsi(i+1,j)-2*tmpsi(i,j)+tmpsi(i-1,j))/(dr^2*R(i,j)) -...
                           (tmpsi(i+1,j)-tmpsi(i-1,j))/(2*dr*R(i,j)^2);
        end
    end
end


% - 固定点磁场
FBR = interp2(R',Z',tmpBR',sep(1,:),sep(2,:));
FBZ = interp2(R',Z',tmpBZ',sep(1,:),sep(2,:));

% - 固定点磁场导数
DFBR = interp2(R',Z',tmpdBR',sep(1,:),sep(2,:));
DFBZ = interp2(R',Z',tmpdBZ',sep(1,:),sep(2,:));

end

