function [BBV,Bpeak] = GradBV(psiout,rh,R,Z)
%   计算真空产磁场分布，及压强峰值位置的磁场大小
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);r1 = R(1,1);
[NR,NZ] = size(R);
ttmpPsi = psiout;
tmpBR = zeros(size(R));tmpBZ = tmpBR;
for i = 1:NR
    for j = 1:NZ
        if ( j>=3 && j<=NZ-2 && i>=3 && i<=NR-2 ) %五点差分
            tmpBR(i,j) = (-1/12*ttmpPsi(i,j+2)+2/3*ttmpPsi(i,j+1)-2/3*ttmpPsi(i,j-1)+...
                          1/12*ttmpPsi(i,j-2))/(dz*R(i,j));
            tmpBZ(i,j) =-(-1/12*ttmpPsi(i+2,j)+2/3*ttmpPsi(i+1,j)-2/3*ttmpPsi(i-1,j)+...
                         1/12*ttmpPsi(i-2,j))/(dr*R(i,j));
        elseif ( j>=2 && j<=NZ-1 && i>=2 && i<=NR-1 )
            tmpBR(i,j)= (ttmpPsi(i,j+1)-ttmpPsi(i,j-1))/(2*dz*R(i,j));
            tmpBZ(i,j)=-(ttmpPsi(i+1,j)-ttmpPsi(i-1,j))/(2*dr*R(i,j)); 
        end
    end
end
BBV = sqrt(tmpBR.^2 + tmpBZ.^2);

Bpeak = BBV( floor((rh - r1)/dr) + 1,find( Z(1,:) == 0) );
end

