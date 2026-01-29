function [BBG,beta,betm] = GradB(psiout,tmpP1,R,Z)
% beta 的最大值在中平面附近检索
mu0 = 4*pi*10^-7;dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);
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
BBG = sqrt(tmpBR.^2 + tmpBZ.^2);
beta = 2*mu0*tmpP1./BBG.^2;
% 计算域边界的NaN转化为0
beta(1,:) = 0;beta(end,:) = 0;beta(:,1) = 0;beta(:,end) = 0;


indz = find( 0 - Z(1,:) <= 0 );% 找距离0最近的点
indzl = indz - floor( 97*2/3 );
indzu = indz + floor( 97*2/3 );
betm = max(max(beta(:,indzl:indzu)));
end

