function [Mo,Ho,Posep] = SIPSNO(tmp_psi,R,Z,Hi,Coils)
% SIPSO Separatrix Inpolygon Points Outer edge (only Separatrix plasma)
% 这里不考虑限制等离子体边界膨胀的限制器；仅包含内侧的线圈的限制器

% update:2024/11/01 更新一个更为稳定的分界面检索算法？
% update:2024/12/22 对稳定性进行更新

% Ver：1.05 仅对分界面和限制器分界线进行区分

% tmp_psi = psi0;
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);dr2 = dr*dr;dz2 = dz*dz;
n  =  10;%收索前20个可疑点

Hir = double(~Hi);  % 内分界面内
NG = length(R(1,:));
RCol = reshape(R,[],1);ZCol=reshape(Z,[],1);

% ---------------------- 外分界面内区域 --------------------------
% 检索区域:上限为L线圈下侧;下限为HC线圈最低点的上侧;右限为HC线圈左侧;左限为D线圈右侧.
%          由输入的Coils线圈位置确定
% 找到半径最小的线圈里最大的半径位置(每个线圈信息的矩阵大小不同)
indR1 = floor((max(max(cell2mat(Coils(2,1))))+0.1 - R(1,1))/dr) + 1;
indR2 = floor((min(min(cell2mat(Coils(3,1))))-0.1 - R(1,1))/dr) + 1;
indZ1 = floor((min(min(cell2mat(Coils(3:end,2))))+0.1 - Z(1,1))/dz) + 1;
indZ2 = floor((max(max(cell2mat(Coils(3:end,2))))-0.1 - Z(1,1))/dz) + 1;
[BBV,~,~] = GradBB(tmp_psi,R,Z);
tgdps = BBV(indR1:indR2,indZ1:indZ2);
tt    = sort(reshape(tgdps,[],1));
tR    = R(indR1:indR2,indZ1:indZ2);
tZ    = Z(indR1:indR2,indZ1:indZ2);
iind  = find(tgdps <= tt(n),n);
rn  = tR(iind);
zn  = tZ(iind);

%%---- 构建的数组 R行不变列变 Z列不变行变
AR = eye(size(rn,1));BR = ones(5,1);KR = kron(AR,BR);
% ---- 位置选取向量
indf  = (1:5:size(KR,1));
inds  = (1:5:size(KR,1))+1;
indt  = (1:5:size(KR,1))+2;
indiv = (1:5:size(KR,1))+3;
indv  = (1:5:size(KR,1))+4;
errs = 1;
toll = 1e-7;
while( errs > 1e-7 )% 五点中心差分
    crr = [rn-2*dr rn-dr rn rn+dr rn+2*dr];
    czz = [zn-2*dz zn-dz zn zn+dz zn+2*dz];
    % - 生产变形矩阵
    CRR = KR*crr;%每三行表示一个点
    CZZ = repmat(reshape(czz',[],1),1,5);
    % 通量插值  行表示y 列表示x mesh与坐标轴一致 ndgrid与(x(行)，y(列))
    VPsi  = interp2(R',Z',tmp_psi',CRR,CZZ,'linear');

    Vdpdr  = (VPsi(indt,1) - 8*VPsi(indt,2) + 8*VPsi(indt,4) -...
        VPsi(indt,5))/(12*dr);
    Vdpdz  = (VPsi(indf,3) - 8*VPsi(inds,3) + 8*VPsi(indiv,3) -...
        VPsi(indv,3))/(12*dz);
    Vdpdr2 = (-VPsi(indt,1) + 16*VPsi(indt,2) -30*VPsi(indt,3) +...
        16*VPsi(indt,4) - VPsi(indt,5))/(12*dr2);
    Vdpdz2 = (-VPsi(indf,3) + 16*VPsi(inds,3) - 30*VPsi(indt,3) +...
        16*VPsi(indiv,3) - VPsi(indv,3))/(12*dz2);
    Vdpdrz =  ((VPsi(indf,1) - 8*VPsi(indf,2) + 8*VPsi(indf,4) - VPsi(indf,5)) -...
        8*(VPsi(inds,1) - 8*VPsi(inds,2) + 8*VPsi(inds,4) - VPsi(inds,5)) +...
        8*(VPsi(indiv,1) - 8*VPsi(indiv,2) + 8*VPsi(indiv,4) - VPsi(indiv,5)) -...
        (VPsi(indv,1) - 8*VPsi(indv,2) + 8*VPsi(indv,4) - VPsi(indv,5)))/(12^2*dr*dz);
    VdetH  = Vdpdr2.*Vdpdz2 - Vdpdrz.^2;

    tmprn = rn + (Vdpdrz.*Vdpdz - Vdpdz2.*Vdpdr)./VdetH;
    tmpzn = zn + (Vdpdrz.*Vdpdr - Vdpdr2.*Vdpdz)./VdetH;
    tmp_VPsi = interp2(R',Z',tmp_psi',tmprn,tmpzn,'linear');
    errs  = abs(max(tmp_VPsi - VPsi(indt,3)));
    rn = tmprn;
    zn = tmpzn;
end
% - 挑选非重复鞍点 ( 里面混杂一个o点 和 两个x点 )
ind_sep = [1;find( abs(diff(tmp_VPsi)) > toll ) + 1];
Vpis_sep = tmp_VPsi(ind_sep);
rn_sep = tmprn(ind_sep);
zn_sep = tmpzn(ind_sep);
re_sep = sign(VdetH(ind_sep));
Posepc = [Vpis_sep';rn_sep';zn_sep';re_sep'];

% - A 区域X点的分界面(选择最小的磁通量)
psi_sepo = min(Posepc(1,find(Posepc(end,:) == -1)));
% ---- 内分界面位置通量
Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo - 1e-5,psi_sepo - 1e-5]);

Posep = Posepc;
% indtt = find(abs(Mo(1,:) - psi_sepo) < abs(min(min(tmp_psi))) &...
%                  Mo(2,:) > 5*max(max(Z)));% 快速定位分界面的每段
indtt = find( Mo(1,:) < 0 );
Ho = 0*R;
for i = 1:length(indtt)
    if ( i ~= length(indtt) )
        tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)) + 1);
    elseif ( i == length(indtt) )
        tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)));
    end
    % 选择区域内包含Hi的区域
    [tIo,~] = inpoly2([RCol ZCol],tMo');
    Io = reshape(tIo,NG,NG);
    indant = Hir.*Io;
    % 除去左上通量 尽量保存右侧通量(这里的作法在于可以有效的分界面划分成有效个小的区域)
    if ( max(max(indant)) == 0 )
        Io = 0*R;
    end
    Ho = Ho + Io;
end

end

