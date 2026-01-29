function [Mo,Mi,Ho,Hi,Posep] = SIPD(tmp_psi,RSOL,RDC,rDC,R,Z)
% SIPD 仅计算D线圈磁通量下的限制器位形的分界面
%      输入内外限值器的位置和磁通量的大小，返回分界面的曲线及Delta函数
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);dr2 = dr*dr;dz2 = dz*dz;
% - 内分界面外区域
psic_tmp1c = interp2(R',Z',tmp_psi',RDC-rDC,0,'linear');
Mi = contourc(R(:,1)',Z(1,:),tmp_psi',[psic_tmp1c,psic_tmp1c]);
NG=length(R(1,:));
RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);
[Mni,~] = inpoly2([RCol ZCol],Mi(:,2:end)');
Ini2 = reshape(Mni,NG,NG);
Hi   = double(~Ini2); % 内分界面外

% ------------------------------------
% - 外分分界面内区域
indzz   = find( Z(1,:) + 50*dz > 0,1) : find( Z(1,:) - 50*dz > 0,1);
psi_lim = min(tmp_psi(floor((RSOL - R(1,1))/dr) + 1,indzz));
Posepl = [psi_lim;RSOL;0;0];
psi_sepo = Posepl(1,1);
% ---- 内分界面位置通量
Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo - 5e-6,psi_sepo - 5e-6]);
[Mno,~] = inpoly2([RCol ZCol],Mo(:,2:end)');
Ino2 = reshape(Mno,NG,NG);
Ho   = double(Ino2); % 内分界面外

% 返回限制器磁通量的LCFS的信息
Posep = [Mi(1,1),Mo(1,1);
         RDC-rDC,RSOL;
         0,0];%第一列内分界面，和第二列外分界面
end

