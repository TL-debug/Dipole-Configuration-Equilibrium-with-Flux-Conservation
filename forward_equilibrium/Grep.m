function [GrP] = Grep(Coils,Sep)
% Sep：等离子体边界点；CRD，CZD D线圈的位置；CRL，CZL L线圈的位置。 
% GREP 计算等离子体边界点与线圈之间的格林函数；
%      包含体积的线圈计算丝电流对边界点的格林函数的平均值。
% 输入边界位置与线圈的位置；输出两者之间的通量值
% NOTE: 平均格林函数的处理大大化简矩阵的维度且是有效的计算固定点和具有匝数线圈之间的方法！！！！（经过测试）
GrP = zeros(size(Coils,1),size(Sep,2));%行表示线圈；列表示等离子体边界位置

for nc = 1:size(Coils,1)
    RC = cell2mat(Coils(nc,1));
    ZC = cell2mat(Coils(nc,2));
    for np = 1:size(Sep,2)
        rpb = Sep(1,np);
        zpb = Sep(2,np);
        
        k2  = 4*RC*rpb./( (RC + rpb).^2 + (ZC - zpb).^2 );
        [MK,ME] = ellipke(k2);
        grp = sqrt(RC.*rpb)./(2*pi*sqrt(k2)).*( (2 - k2 ).*MK - 2*ME );

        GrP(nc,np) = mean(mean(grp));
%         GrP(nc,np) = sum(sum(grp));
    end
end
GrP = -GrP';
end

% %% - 在靠近D线圈位置的边界使用格林函数计算的磁通量具有0.05%的相对误差较大，
% % 绘制对应点位置
% plot(Mo(1,2:end),Mo(2,2:end),'.m');hold on;
% plot(Sep(1,:),Sep(2,:),'.k','markersize',25);axis([0 4 -1.5 2.0])
% for k = 1:length(Sep)
%    text(Sep(1,k)+0.1,Sep(2,k)+0.1,[num2str(k)],'fontsize',14) ;
% end
% set(gca,'FontName','Times','FontSize',18);
% set(gca,'TickLabelInterpreter','latex');
% xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",18);
% ylabel('$Z\ (\rm m)$','Interpreter','latex',"FontSize",18);
% 
% axes('position',[0.37 0.36 0.4 0.30]);
% intpsiD = interp2(R',Z',PSICD',Sep(1,:),Sep(2,:));
% intpsiL = interp2(R',Z',PSICL',Sep(1,:),Sep(2,:));
% GrepsiD = GrP(2,:)*CJ(2)*mu0;
% GrepsiL = GrP(1,:)*CJ(1)*mu0;
% % 两种算法的误差对比
% semilogy(abs(intpsiD - GrepsiD)./abs(intpsiD),'-b','linewidth',2);hold on;
% semilogy(abs(intpsiL - GrepsiL)./abs(intpsiL),'--r','linewidth',2);axis tight;
% grid on;ylabel error;xlabel point;