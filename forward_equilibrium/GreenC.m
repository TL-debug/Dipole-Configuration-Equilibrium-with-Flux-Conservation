function [GrC,indjb,nb] = GreenC(Coils,rr,zz)
% GreenC 计算格林函数在计算域边界与线圈之间的格林函数
%      返回的一个二维矩阵行边界点，列计算域网格
%      这里可以使用将计算域的电流变形为一列，使用矩阵乘一次性求解所有边界点磁通量，
%      再根据位置指引向量对边界进行赋值
% 使用规则：tpsibc = -mu0*GrE*reshape(RHSC,[],1);
%    Note：每一行表示计算域边界位置；列表示线圈位置。默认计算域不包含R=0。

NR = length(rr);NZ = length(zz);
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


[R,Z] = ndgrid(rr,zz);% 先计算后转为向量
GrC = zeros(length(jbr),size(Coils,1));
tR = cell2mat(Coils(:,1))';% 线圈位置一行
tZ = cell2mat(Coils(:,2))';%

for k1 = 1:nb
     k2 = 4*tR*R(jbr(k1),jbz(k1))./( (R(jbr(k1),jbz(k1)) + tR).^2 +...
                                     (Z(jbr(k1),jbz(k1)) - tZ).^2 );
    [MK,ME] = ellipke(k2);

    GrC(k1,:) = sqrt(tR.*R(jbr(k1),jbz(k1)))./(2*pi*sqrt(k2)).*...
                                         ( (2 - k2).*MK - 2*ME );
end

indjb = [jbr;jbz]; % jbr和jbz都是行向量
end

