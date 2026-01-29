function [GrE,indjb,nb] = Gree(rr,zz)
% Gree 计算边界点与计算域网格的格林函数，返回的一个二维矩阵行边界点，列计算域网格
%      这里可以使用将计算域的电流变形为一列，使用矩阵乘一次性求解所有边界点磁通量，
%      再根据位置指引向量对边界进行赋值

NR = length(rr);NZ = length(zz);
if ( rr(1) == 0 )% 包含R=0
    % 三边界
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
elseif ( rr(1) ~= 0 )
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
end

[R,Z] = ndgrid(rr,zz);% 先计算后转为向量
tR = R(2:end-1,2:end-1);
tZ = Z(2:end-1,2:end-1);
GrE = zeros(nb,length(reshape(R,[],1)));
TGre = zeros(size(R));
for k1 = 1:nb
    k2 = 4*tR*R(jbr(k1),jbz(k1))./( (R(jbr(k1),jbz(k1)) + tR).^2 +...
                                    (Z(jbr(k1),jbz(k1)) - tZ).^2 );
    [MK,ME] = ellipke(k2);
    tGre = sqrt(tR.*R(jbr(k1),jbz(k1)))./(2*pi*sqrt(k2)).*...
                                         ( (2 - k2).*MK - 2*ME );
    TGre(2:end-1,2:end-1) = tGre;
    GrE(k1,:) = reshape(TGre,[],1)';
end

indjb = [jbr;jbz]; % jbr和jbz都是行向量
end


% %% 对GRE1和GRE1进行测试
% psibc1 = -mu0*GrE1*reshape(tmpJ1,[],1)*ds;
% 
% psibc2 = 0*psibc1;
% for j = 1:nb
%     psibc2(j) = -mu0*ds*reshape(GrE(:,:,j),[],1)'*reshape(tmpJ1,[],1);
% end
% plot(psibc1,'-b','linewidth',4);hold on;
% plot(psibc2,'--r','linewidth',4);