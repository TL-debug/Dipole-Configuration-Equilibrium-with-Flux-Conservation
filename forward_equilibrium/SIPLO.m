function [Mo,Ho,Posep] = SIPLO(tmp_psi,RSOL,R,Z,Hi)
% SIPL 仅在一开始计算，后续则使用其他分界面子程序
% a.分界面
% b.限制器
% - 测试用
%tmp_psi = psi0;
% Ver：1.05 仅对分界面和限制器分界线进行区分
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);dr2 = dr*dr;dz2 = dz*dz;
Hir  = double(~Hi);  % 内分界面内
NG=length(R(1,:));
RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);

% - 外分界面 ( 自动检索 )
n  =  10;     % 收索前20个可疑点
% 在指定区域内检索
indr = find( R(:,1) >= 0.9 & R(:,1) <= 1.6 );
indz = find( Z(1,:) >= 1.3 & Z(1,:) <= 2.0 );
tR = R(indr,indz);tZ = Z(indr,indz);
[BBV,~,~] = GradBB(tmp_psi(indr,indz),tR,tZ);
tgdps = BBV(2:end-1,2:end-1);
%mesh(tR(2:end-1,2:end-1),tZ(2:end-1,2:end-1),tgdps)
tt    = sort(reshape(tgdps,[],1));
iind  = find(tgdps <= tt(n),n);
vr    = tR(2:end-1,2:end-1);vz = tZ(2:end-1,2:end-1);
rn    = vr(iind);zn = vz(iind);
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

% - A 区域X点的分界面
psi_sepo = min(Posepc(1,find(Posepc(end,:) == -1)));
Posep = Posepc;
% ---- 内分界面位置通量
Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo,psi_sepo] - 5e-6);

if ( max(Mo(1,:)) > RSOL )
    indzz   = find( Z(1,:) + 50*dz > 0,1) : find( Z(1,:) - 50*dz > 0,1);
    psi_lim = min(tmp_psi(floor((RSOL - R(1,1))/dr) + 1,indzz));
    Posepl = [psi_lim;RSOL;0;0];
    Posep = Posepl;
    psi_sepo = Posepl(1,1);
    % ---- 内分界面位置通量
    Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo,psi_sepo]);
end


indtt = find(abs(Mo(1,:) - psi_sepo) < abs(min(min(tmp_psi))) &...
                 Mo(2,:) > 5*max(max(Z)));% 快速定位分界面的每段
Ho = 0*R;
for i = 1:length(indtt)
    if ( i ~= length(indtt) )
        tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)) + 1);
    elseif ( i == length(indtt) )
        tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)));
    end
    [tIo,~] = inpoly2([RCol ZCol],tMo');
    Io = reshape(tIo,NG,NG);
    indant = Hir.*Io;
    % 除去左上通量 尽量保存右侧通量
    if ( max(max(indant)) == 0 )
        Io = 0*R;
    end
    Ho = Ho + Io;
end




% % Ver：1.1
% dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);dr2 = dr*dr;dz2 = dz*dz;
% n  =  20;%收索前20个可疑点
% 
% Hir  = double(~Hi);  % 内分界面内
% NG=length(R(1,:));
% RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);
% % - 外分界面
% [BBV,~,~] = GradBB(tmp_psi,R,Z);
% tgdps = BBV(5:end-4,5:end-4);
% tt    = sort(reshape(tgdps,[],1));
% tR    = R(5:end-4,5:end-4);
% tZ    = Z(5:end-4,5:end-4);
% iind = find(tgdps <= tt(n),n);
% vr  = tR(iind);vz  = tZ(iind);
% vr1 = vr;vz1 = vz;vr2 = vr;vz2 = vz;
% % - 选择A区域可疑点（随机5个）
% indx = find( vr1 < 0.9 | vr1 > 1.9 );vr1( indx ) = [];vz1( indx ) = [];
% indy = find( vz1 < 1   | vz1 > 2   );vz1( indy ) = [];vr1( indy ) = [];
% ind15 = floor( rand(1,5)*(length(vr1) - 1)) + 1;
% % - 选择B区域可疑点（随机5个）
% indx = find( vr2 < 1    | vr2 > RSOL);vr2( indx ) = [];vz2( indx ) = [];
% indy = find( vz2 < -0.2 | vz2 > 0.5  );vz2( indy ) = [];vr2( indy ) = [];
% % ind25 = floor( rand(1,5)*(length(vr2) - 1)) + 1;
% 
% rn = vr1(ind15);zn = vz1(ind15);
% %%---- 构建的数组 R行不变列变 Z列不变行变
% AR = eye(size(rn,1));BR = ones(5,1);KR = kron(AR,BR);
% % ---- 位置选取向量
% indf  = (1:5:size(KR,1));
% inds  = (1:5:size(KR,1))+1;
% indt  = (1:5:size(KR,1))+2;
% indiv = (1:5:size(KR,1))+3;
% indv  = (1:5:size(KR,1))+4;
% errs = 1;
% toll = 1e-7;
% while( errs > 1e-7 )% 五点中心差分
%     crr = [rn-2*dr rn-dr rn rn+dr rn+2*dr];
%     czz = [zn-2*dz zn-dz zn zn+dz zn+2*dz];
%     % - 生产变形矩阵
%     CRR = KR*crr;%每三行表示一个点
%     CZZ = repmat(reshape(czz',[],1),1,5);
%     % 通量插值  行表示y 列表示x mesh与坐标轴一致 ndgrid与(x(行)，y(列))
%     VPsi  = interp2(R',Z',tmp_psi',CRR,CZZ,'linear');
%  
%     Vdpdr  = (VPsi(indt,1) - 8*VPsi(indt,2) + 8*VPsi(indt,4) -...
%               VPsi(indt,5))/(12*dr);
%     Vdpdz  = (VPsi(indf,3) - 8*VPsi(inds,3) + 8*VPsi(indiv,3) -...
%               VPsi(indv,3))/(12*dz);
%     Vdpdr2 = (-VPsi(indt,1) + 16*VPsi(indt,2) -30*VPsi(indt,3) +...
%               16*VPsi(indt,4) - VPsi(indt,5))/(12*dr2);
%     Vdpdz2 = (-VPsi(indf,3) + 16*VPsi(inds,3) - 30*VPsi(indt,3) +...
%               16*VPsi(indiv,3) - VPsi(indv,3))/(12*dz2);
%     Vdpdrz =  ((VPsi(indf,1) - 8*VPsi(indf,2) + 8*VPsi(indf,4) - VPsi(indf,5)) -...
%              8*(VPsi(inds,1) - 8*VPsi(inds,2) + 8*VPsi(inds,4) - VPsi(inds,5)) +...
%              8*(VPsi(indiv,1) - 8*VPsi(indiv,2) + 8*VPsi(indiv,4) - VPsi(indiv,5)) -...
%             (VPsi(indv,1) - 8*VPsi(indv,2) + 8*VPsi(indv,4) - VPsi(indv,5)))/(12^2*dr*dz);
%     VdetH  = Vdpdr2.*Vdpdz2 - Vdpdrz.^2;
%     
%     tmprn = rn + (Vdpdrz.*Vdpdz - Vdpdz2.*Vdpdr)./VdetH;
%     tmpzn = zn + (Vdpdrz.*Vdpdr - Vdpdr2.*Vdpdz)./VdetH;
%     tmp_VPsi = interp2(R',Z',tmp_psi',tmprn,tmpzn,'linear');
%     errs  = abs(max(tmp_VPsi - VPsi(indt,3)));
%     rn = tmprn;
%     zn = tmpzn;
% end
% % - 挑选非重复鞍点 ( 里面混杂一个o点 和 两个x点 )
% ind_sep = [1;find( abs(diff(tmp_VPsi)) > toll ) + 1];
% Vpis_sep = tmp_VPsi(ind_sep);
% rn_sep = tmprn(ind_sep);
% zn_sep = tmpzn(ind_sep);
% re_sep = sign(VdetH(ind_sep));
% Posepc = [Vpis_sep';rn_sep';zn_sep';re_sep'];
% 
% % - A 区域X点的分界面
% psi_sepo = min(Posepc(1,find(Posepc(end,:) == -1)));
% % ---- 内分界面位置通量
% Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo - 5e-6,psi_sepo - 5e-6]);
% 
% if ( max(Mo(1,:)) > RSOL )
% 	% - 限制器位置的通量
%     indzz   = find( Z(1,:) + 80*dz > 0,1) : find( Z(1,:) - 80*dz > 0,1);
%     psi_lim = min(tmp_psi(floor((RSOL - R(1,1))/dr) + 1,indzz));
%     Posepl = [psi_lim;RSOL;0;0];
%     
%     if ( isenpty( vr2 ) == 1 )
%         % - A 区域X点的分界面
%         psi_sepo = Posepl(1,1);
%         % ---- 内分界面位置通量
%         Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo - 5e-6,psi_sepo - 5e-6]);
%         
%     else % 区域B与限制器分界面判断
%         ind25 = floor( rand(1,5)*(length(vr2) - 1)) + 1;
%         ind25 = [1,ind25( find( abs(diff(ind25)) > 0.01 ) +1)];
%         % - 计算区域B可疑点
%         vr2 = vr2(ind25);vz2 = vz2(ind25);% 非重复的指标
%         Posepp = zeros(4,length(vr1));
%         
%         for k = 1:length(vr2)
%             rn = vr2(k);zn = vz2(k);
%             errs = 1;
%             while ( errs > toll )%此处第五次迭代数值误差小于10^-14
%                 crr=[rn-dr rn rn+dr];
%                 czz=[zn-dz zn zn+dz];
%                 [CRR,CZZ]=ndgrid(crr,czz);
%                 Psic_tmp=interp2(R',Z',tmp_psi',CRR',CZZ','line');
%                 dpdr=[0 1 0]*Psic_tmp*[-1 0 1]'/(2*dr);
%                 dpdz=[-1 0 1]*Psic_tmp*[0 1 0]'/(2*dz);
%                 dpdr2=[0 1 0]*Psic_tmp*[1 -2 1]'/dr2;
%                 dpdz2=[1 -2 1]*Psic_tmp*[0 1 0]'/dz2;
%                 dpdrdz=[-1 0 1]*Psic_tmp*[-1 0 1]'/(4*dr*dz);
%                 detH=dpdr2*dpdz2-dpdrdz^2;
%                 rn1=rn+(dpdrdz*dpdz-dpdz2*dpdr)/(detH);
%                 zn1=zn+(dpdrdz*dpdr-dpdr2*dpdz)/(detH);
%                 %新位置的插值
%                 psic_tmp1=interp2(R',Z',tmp_psi',rn1,zn1,'line');
%                 errs = abs(psic_tmp1-Psic_tmp(2,2));
%                 rn=rn1;
%                 zn=zn1;
%             end
%             Posepp(1,k) = psic_tmp1;
%             Posepp(2,k) = rn;
%             Posepp(3,k) = zn;
%             Posepp(4,k) = sign(detH);
%         end
%         % 删除相同值及无效值（B区域有可疑点，迭代后）
%         % - 删除无效值
%         Posepp(:,find( Posepp(2,:) < 0.5 | Posepp(2,:) > RSOL )) = [];
%         Posepp(:,find( Posepp(3,:) < -1 | Posepp(3,:) > 1 )) = [];
%         Posepp(:,isnan(Posepp(1,:))) = [];
%         % - 删除相同值
%         if ( isempty(Posepp) == 0 )
%             ind_sep = [1;find( abs(diff(Posepp(1,:))) > toll )+1];
%             Posepp = Posepp(:,ind_sep);
%         end
%         % - 挑选分界面区域
%         if ( isempty(Posepp) == 1 || isempty(Posepp(4,:) == -1) == 1 ) % 删除o点
%             psi_sepo = Posepl(1,1);
%         elseif ( isempty(Posepp(4,:) == -1) == 0 ) % x点 选择通量值小的为边界
%             % 找Posepp中最小的
%             psi_sepo = min(min(Posepp(1,find(Posepp(end,:) == -1))),Posepl(1,1));
%         end
%         Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo - 5e-6,psi_sepo - 5e-6]);
%     end
% 
% end
% 
% indtt = find(abs(Mo(1,:) - psi_sepo) < abs(min(min(tmp_psi))) &...
%     Mo(2,:) > 5*max(max(Z)));% 快速定位分界面的每段
% Ho = 0*R;
% for i = 1:length(indtt)
%     if ( i ~= length(indtt) )
%         tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)) + 1);
%     elseif ( i == length(indtt) )
%         tMo = Mo( :,indtt(i) + 1 : indtt(i) + Mo(2,indtt(i)));
%     end
%     
%     [tIo,~] = inpoly2([RCol ZCol],tMo');
%     Io = reshape(tIo,NG,NG);
%     indant = Hir.*Io;
%     % 除去左上通量 尽量保存右侧通量
%     if ( max(max(indant)) == 0 )
%         Io = 0*R;
%     end
%     Ho = Ho + Io;
% end

end
