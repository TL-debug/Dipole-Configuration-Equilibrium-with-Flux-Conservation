function [Mo,Ho,Posep] = SIPDO(tmp_psi,RSOL,R,Z,Hi)
%SIPDO 仅计算D线圈磁通量下的限制器位形的分界面，仅更新outer边界的形状和数值
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);dr2 = dr*dr;dz2 = dz*dz;
Hir  = double(~Hi);  % 内分界面内
NG=length(R(1,:));
RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);

indzz   = find( Z(1,:) + 50*dz > 0,1) : find( Z(1,:) - 50*dz > 0,1);
psi_lim = min(tmp_psi(floor((RSOL - R(1,1))/dr) + 1,indzz));
Posepl = [psi_lim;RSOL;0];
Posep = Posepl;
psi_sepo = Posepl(1,1);
% ---- 内分界面位置通量
Mo = contourc(R(:,1)',Z(1,:),tmp_psi',[psi_sepo,psi_sepo]);

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

end
% % ---- 测试等离子体区域识别代码
% [Mo,Mi,Ho,Hi,Posep] = SIPD(PSIC,RSOL,RDC,rDC,R,Z);
% 
% plot(Mi(1,2:end),Mi(2,2:end),'-m','linewidth',2);hold on;
% plot(Mo(1,2:end),Mo(2,2:end),'-m','linewidth',2);
% axis([r1 r2 z1 z2]);
% 
% [Mo,Ho,Posep] = SIPDO(PSIC,RSOL,R,Z,Hi);
% plot(Mi(1,2:end),Mi(2,2:end),'--g','linewidth',2);hold on;
% plot(Mo(1,2:end),Mo(2,2:end),'--g','linewidth',2);