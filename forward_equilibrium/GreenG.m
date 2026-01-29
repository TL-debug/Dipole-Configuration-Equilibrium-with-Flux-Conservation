function [GrE,jbr,jbz,nb] = GreenG(rr,zz)
% GreenG 计算域边界对计算域网格的格林函数矩阵
% 功能：1.R = 0时，格林函数仅计算3个边界
%       2.R ~ 0时，格林函数计算域中四个边界
NR = length(rr);NZ = length(zz);
if ( rr(1) == 0 )% 包含R=0

    nb = NZ + 2*NR - 4;
    % ----------- 定位边界 ----------- %
    % 下边界
    jbr1 = 2:NR;
    jbz1 = ones(1,length(jbr1));
    % 右边界
    jbz2 = 2:NZ-1;
    jbr2 = NR*ones(1,length(jbz2));
    % 上边界
    jbr3 = jbr1;
    jbz3 = NZ*ones(1,length(jbr3));
    % ------------- 边界位置 ------------ %
    jbr = [jbr1 jbr2 jbr3];
    jbz = [jbz1 jbz2 jbz3];
    
    GrE = zeros(NR,NZ,nb);
    [R,Z] = ndgrid(rr,zz);
    for i = 1:nb
        k2 = 4*R(2:end-1,2:end-1)*R(jbr(i),jbz(i))./((R(2:end-1,2:end-1)+...
            R(jbr(i),jbz(i))).^2+(Z(2:end-1,2:end-1)-Z(jbr(i),jbz(i))).^2);
        [MK,ME] = ellipke(k2);
        GrE(2:end-1,2:end-1,i) = sqrt(R(2:end-1,2:end-1)*R(jbr(i),jbz(i)))./...
            (2*pi*sqrt(k2)).*((2-k2).*MK-2*ME);
    end
    
elseif ( rr(1) ~= 0 ) % 不包含R = 0
    % 四个边界

    nb = 2*NR + 2*NZ -4;
    % 将边界ind进行逆时针排列；（逆时针方向，曲线积分的正方向）
    % 下边界
    jbr1 = 2:NR-1;
    jbz1 = ones(1,length(jbr1));
    % 右边界
    jbz2 = 1:NZ;
    jbr2 = NR*ones(1,length(jbz2));
    % 上边界
    jbr3 = flip(jbr1);
    jbz3 = NZ*ones(1,length(jbr3));
    % 左边界
    jbz4 = flip(jbz2);
    jbr4 = ones(1,length(jbz4));
    % 边界的indx
    jbr = [jbr1 jbr2 jbr3 jbr4];
    jbz = [jbz1 jbz2 jbz3 jbz4];
    
    GrE = zeros(NR,NZ,nb);
    [R,Z] = ndgrid(rr,zz);
    for i = 1:nb
        k2 = 4*R(2:end-1,2:end-1)*R(jbr(i),jbz(i))./((R(2:end-1,2:end-1)+...
            R(jbr(i),jbz(i))).^2+(Z(2:end-1,2:end-1)-Z(jbr(i),jbz(i))).^2);
        [MK,ME] = ellipke(k2);
        GrE(2:end-1,2:end-1,i) = sqrt(R(2:end-1,2:end-1)*R(jbr(i),jbz(i)))./...
            (2*pi*sqrt(k2)).*((2-k2).*MK-2*ME);
    end
    
end
end