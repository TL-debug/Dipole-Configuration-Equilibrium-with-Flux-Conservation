function [MBT] = minSB(tmp_psi1,BBG,Hi,Ho,R,Z)
% MINBB 输入通量和磁场强度 返回各点沿磁场线的最小值用于计算压强各项异性的平衡
%   SB is A Symmetric configuration of magnetic field B
% dr = R(2,1) - R(1,1);
[NR,NZ] = size(R);

% ---- 沿磁场线检索
n = 1;
% 采样区域（在中平面上）
indz = find( Z(1,:) == 0 );
% 内边界之外
indr1 = max(find( Ho(:,indz) == 1 ) );
% 外边界之内
indr2 = max(find( Hi(:,indz) == 0 )) + 1;

% 均匀采样点前1/4区域取采样点密集，后3/4采样点稀疏
sd = 5;
tindr = floor(abs(indr2 - indr1)/sd) + indr2;
indr = [indr2:n:tindr,tindr+1:1*n:indr1];
% 所有采样点的磁场线位置向量
spp = tmp_psi1(indr,indz);
M = contourc(R(:,1)',Z(1,:),tmp_psi1',spp);% 是否可以返回指定区域磁场线，并找最小值，可大大减少数据量
tmpM = M;
% 目的将M中位置相关信息去掉，并重新存储位置信息到另一个矩阵，M剩下纯R和Z点信息
% 1先将原始位置信息挑选出来
tmpind = find( M(1,:) < 0  & M(2,:) > 10 );
indd(:,1) = tmpind'+1 ;% 第一列为起始位置；第二列为终点位置
indd(:,2) = [tmpind(2:end)'-1;length(M(1,:))];
indnn = 1:length(spp);
indnn = repmat(indnn',1,2);
indd = indd - indnn; % 采样点对应的起始和终止位置的坐标
tmpM(:,tmpind) = [];% 所有采样点磁场线的位置

VB = interp2(R',Z',BBG',tmpM(1,:),tmpM(2,:),'linear');%所有采样点的对应的磁场强度
% 选取磁场线区域间的磁场最小值并找到最小值位置
MBB  = zeros(size(spp));% 沿磁场线最小值
for j4 = 1:size(indd,1)
    tmpBB = VB(indd(j4,1):indd(j4,2));
%     tmprr = tmpM(1,indd(j4,1):indd(j4,2));
%     tmpzz = tmpM(2,indd(j4,1):indd(j4,2));
    % 找到最小值
    [minB,~] = min(tmpBB);
    % 返回反之
    MBB(j4)    = minB; 
%     MBRZ(:,j4) =[tmprr(tmpindd);tmpzz(tmpindd)];
end


MBT = reshape(0*R,[],1);
PVT = reshape(tmp_psi1,[],1);
MBT = interp1(spp,MBB,PVT,'nearest');
MBT = reshape(MBT,NR,NZ);
% MBT(MBT == 0) = NaN;
% % - 由于对称位形关于z=0，假设最小值在中平面上
% n = 1;
% % 采样区域（在中平面上）
% indz = find( Z(1,:) == 0 );
% % 内边界之外
% indr1 = max(find( Ho(:,indz) == 1 ) );
% % 外边界之内
% indr2 = max(find( Hi(:,indz) == 0 )) + 1;
% tic
% % 均匀采样点
% indr = indr2:n:indr1;
% mbr  = R(indr,indz);
% % 所有采样点的磁场线位置向量
% spp = tmp_psi1(indr,indz);
% 
% 
% 
% 
% VB  = 0*spp;% 采样点沿磁场线的磁场最小值
% 
% VB(2:end-1) = (spp(3:end) - spp(1:end-2))./(2*dr*mbr(2:end-1));
% VB(1)       = (spp(2) - spp(1))./(dr*mbr(1));
% VB(end)     = (spp(end) - spp(end-1))./(dr*mbr(end));
% 
% 
% % - 对通量进行插值计算沿磁场线的磁场最小值
% BMT = reshape(0*R,[],1);
% PVT = reshape(tmp_psi1,[],1);
% MBT = interp1(spp,VB,PVT,'nearest');
% MBT = reshape(MBT,NR,NZ);
end

