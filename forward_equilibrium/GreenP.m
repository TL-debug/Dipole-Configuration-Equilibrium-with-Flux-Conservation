% function [Grp] = GreenP(R,Z,Sep)
% % ---- 假设所有源都位于计算域边界内
% % Ver1.0:GREENP 计算等离子体边界点对计算域网格的格林函数
% % 输入 计算域网格 与 等离子体边界点的位置 ； 输出 两者之间的格林函
% Grp = zeros(length(R(:,1)),length(Z(1,:)),length(Sep(1,:)));
% tR = R(2:end-1,2:end-1);% 这样处理为了方便处理计算域的R=0的情况
% tZ = Z(2:end-1,2:end-1);
% for k = 1:length(Sep(1,:))
%     k2 = 4*tR*Sep(1,k)./( (tR + Sep(1,k)).^2 + (tZ - Sep(2,k)).^2 );
%     [MK,ME] = ellipke(k2); 
%     tmp_Grp = sqrt(tR.*Sep(1,k))./(2*pi*sqrt(k2)).*( (2 - k2).*MK - 2*ME );
%     if ( max ( tmp_Grp == inf ) == 1 )
%         ind_Grp = find( tmp_Grp == inf );
%         Grp_int = interp2(tR',tZ',tmp_Grp',tR(ind_Grp),tZ(ind_Grp));
%         tmp_Grp( ind_Grp ) = Grp_int;
%     end
%     Grp(2:end-1,2:end-1,k) = tmp_Grp;
% end

function [Grp,Mcs,CCoils] = GreenP(Coils,Sep)
% ---- 假设所有源都位于计算域边界内,且线圈最小半径不能为0
% Ver2.0:GREENP 计算等离子体边界点对线圈丝电流之间的格林函数
% 输入 线圈的信息 与 等离子体边界点的位置;note:这里不分DL线圈的区别，但是线圈电流要对应正确
% 输出 矩阵Grp行表示线圈，列表示边界点，一个将丝电流组合成线圈电流的矩阵，一个初始丝电流向量

% 返回电流丝的R和Z方向位置向量，并返回对应的电流丝电流向量
CR = []; CZ = []; CI = []; FCN = [];
for k1 = 1:size(Coils,1)
    tr = reshape(cell2mat(Coils(k1,1)),[],1);
    tz = reshape(cell2mat(Coils(k1,2)),[],1);
    ti = cell2mat(Coils(k1,3))/length(tz)*ones(size(tz));

    CR = [CR;tr];
    CZ = [CZ;tz];
    CI = [CI;ti];
    FCN = [FCN,length(tr)];
end
% 合成矩阵
Mcs = [];
for k1 = 1:length(FCN)
    tmcs = ones(1,FCN(k1));
    Mcs  = blkdiag(Mcs,tmcs);
end

% 线圈丝电流位置和电流信息
CCoils = [CR,CZ,CI];

% 丝电流和位置间的格林函数
Grp = zeros(length(CR),length(Sep(1,:)));
for k1 = 1:length(Sep(1,:))
    k2 = 4*Sep(1,k1)*CR./( (CR + Sep(1,k1)).^2 + (CZ - Sep(2,k1)).^2 );
    [MK,ME] = ellipke(k2);
    Grp(:,k1) = 0.5/pi*sqrt(Sep(1,k1)*CR./k2).*( (2 - k2).*MK - 2*ME );
end

end

% ---- 使用该子程序的写法
% [GrP,Mcs,CCoils] = GreenP(Coils,Sep);
% % 计算固定点上的磁通量（测试）
% psifpc = zeros(size(Sep(1,:)));
% for k1 = 1:length(Sep(1,:))
%     psifpc(k1) = -mu0*GrP(:,k1)'*CCoils(:,3);
% end




%% ---- 参考互感算法
% for k1 = 1:NCD
%     % k2不计算R=0的点，
%     k2 = 4*VRD(k1)*tRR./( (tRR + VRD(k1)).^2 + (tZZ - VZD(k1)).^2 );
%     [MK,ME] = ellipke(k2);
%     tmp_MuT = mu0*sqrt(tRR.*VRD(k1)./k2).*( (2 - k2).*MK - 2*ME );
%     if ( max ( tmp_MuT == Inf ) == 1 )
%         % - 插值计算无穷点
%         ind_Minf = find( tmp_MuT == inf );
%         mut_int  = interp2(tRR',tZZ',tmp_MuT',tRR(ind_Minf),tZZ(ind_Minf));
%         tmp_MuT ( ind_Minf ) = mut_int;
%     end
%     MuT(2:end-1,2:end-1,k1) = tmp_MuT;
% end
%% ---- 测试格林函数计算等离子体通量在等离子体边界点数值
% 检查边界点磁通量是否计算正确 test ok
% plot3(Sep(1,:),Sep(2,:),psicp,'or','linewidth',2,'markersize',10);hold on;
% plot(Sep(1,:),Sep(2,:),'.g','markersize',20)
% plot(Mo(1,2:end),Mo(2,2:end),'--m');
% mesh(R,Z,PSIC);
% axis([r1 r2 z1 z2]);

%% ---- 对比平均格林函数和丝电流计算固定点磁通量的对比
% 
% % - 猜测初始值
% clear;clc;close all;
% % 构建计算域网格
% nI = 8;NG = 2^nI + 1;nlay = nI - 1;%多重网格层数
% NR = NG;NZ = NG;c0 = 2e-7;r1 = 0.0;r2 = 4.0;z1 = -1.5;z2 = 2.5;mu0 = 4*pi*1e-7;
% rr=linspace(r1,r2,NR);zz=linspace(z1,z2,NZ);dr=rr(2)-rr(1);dz=zz(2)-zz(1);
% dr2=dr*dr;dz2=dz*dz;ds=dr*dz;
% %大线圈几何参数(中心位置和拆解个数）
% CR = [0.85 0.53];CZ = [2.1 0];LNR = 60;LNZ = 30;DNR = 40;DNZ = 40;
% %悬浮大线圈的参数
% LW = 0.3;LH = 0.15;
% %偶极场大线圈的参数
% DW = 0.16;DH = 0.16;
% % 线圈电流 ( Current of Coils )
% CJ = [4.7e5 4.77e6];
% %偶极场线圈
% VCRD = linspace(CR(2) - DW/2,CR(2) + DW/2,DNR);
% VCZD = linspace(CZ(2) - DH/2,CZ(2) + DH/2,DNZ);
% %悬浮线圈
% VCRL = linspace(CR(1) - LW/2,CR(1) + LW/2,LNR);
% VCZL = linspace(CZ(1) - LH/2,CZ(1) + LH/2,LNZ);
% [CRD,CZD] = ndgrid(VCRD,VCZD);
% [CRL,CZL] = ndgrid(VCRL,VCZL);
% % 写为线圈元的L和D线圈位置R和Z的位置
% Coils = [{CRL},{CZL},{CJ(1)};  % L_Coil
%          {CRD},{CZD},{CJ(2)}]; % D_Coil   
% %------生成空间中的磁场和通量子程序
% [R,Z]=ndgrid(rr,zz);
% [PSICD] = MCVD(rr,zz,CRD,CZD,CJ(2));
% [PSICL] = MCVD(rr,zz,CRL,CZL,CJ(1));
% PSIC = PSICD + PSICL;
% % ----- 生成插值和弛豫矩阵（R和Z）
% [MRI] = MatRI(R,Z,nlay);        % 返回矩阵中包括插值和弛豫矩阵及每层网格
% % ---------- 生成格林函数矩阵GrE ------------- %
% [GrE,nbr,nbz,nb] = GreenG(rr,zz);
% psibct  = zeros(1,nb); % total magnetic flux of the computational domain
% psibcd = zeros(1,nb);
% psibcl = zeros(1,nb);
% for j = 1:nb
% 	psibct(j)  = PSIC(nbr(j),nbz(j));
%     psibcd(j) = PSICD(nbr(j),nbz(j));
%     psibcl(j) = PSICL(nbr(j),nbz(j));
% end
% % ------ p2gv 配置求解区域椭圆算符的RHS项 ------ %
% [RHSCt] = p2gv(Coils,rr,zz); 
% RHSCt = mu0*R.*RHSCt/ds; % total source contribution (RHS)
% % ----- D线圈自感，D与P互感，D与L的互感
% [mDL] = FMI(CRD,CZD,CRL,CZL);  % 丝电流方法互感
% [InD] = FSI(CRD,CZD);          % 丝电流方法自感
% [ncd,MuT] = GMI(R,Z,CRD,CZD);  % D线圈电流丝域网格间的互感
% % ------ 读取等离子体边界信息
% load('sep_points.mat');% 读取等离子体分界面的点的位置信息
% % 计算线圈丝电流对固定点的格林函数, Grp 行为丝电流格林函数，列为固定点。
% % Grp 和 Coils将线圈按照Coil中的顺序展开的为一列向量的形式。
% [GrP1] = Grep(Coils,Sep);
% [GrP,Mcs,CCoils] = GreenP(Coils,Sep);
% IDLP0 = IDLP;
% % 计算固定点上的磁通量（测试）
% psifp = zeros(size(Sep(1,:)));
% for k1 = 1:length(Sep(1,:))
%     psifp(k1) = -mu0*GrP(:,k1)'*CCoils(:,3);
% end
% 
% subplot(2,2,1)
% tpsifp = interp2(R',Z',PSIC',Sep(1,:),Sep(2,:),'cubic');
% GrP2(1,:) = mu0*GrP1(1,:)*CJ(1);
% GrP2(2,:) = mu0*GrP1(2,:)*CJ(2);
% psifp1 = sum(GrP2);
% plot(psifp1,'or','markersize',10,'linewidth',2);hold on;
% plot(tpsifp,'.k','markersize',20);
% box on;grid on;
% indmaxerr = find(abs( (tpsifp - psifp)./psifp ) == max(abs( (tpsifp - psifp)./psifp )));
% plot(indmaxerr,psifp(indmaxerr),'pg','markersize',20,'linewidth',1);
% title({'Average Green Function';['Max Relative Error = ',num2str(max(abs( (tpsifp - psifp)./psifp )))]})
% xlabel('$points\ number$','Interpreter','latex',"FontSize",16);
% ylabel('$\psi\ (\rm wb)$','Interpreter','latex',"FontSize",16);
% set(gca,'FontName','Times','FontSize',16);
% 
% subplot(2,2,2)
% plot(psifp,'or','markersize',10,'linewidth',2);hold on;
% plot(tpsifp,'.k','markersize',20);
% 
% indmaxerr = find(abs( (tpsifp - psifp)./psifp ) == max(abs( (tpsifp - psifp)./psifp )));
% plot(indmaxerr,psifp(indmaxerr),'pg','markersize',20,'linewidth',1);
% box on;grid on;
% title({'Current Filament method';['Max Relative Error = ',num2str(max(abs( (tpsifp - psifp)./psifp )))]})
% xlabel('$points\ number$','Interpreter','latex',"FontSize",16);
% ylabel('$\psi\ (\rm wb)$','Interpreter','latex',"FontSize",16);
% set(gca,'FontName','Times','FontSize',16);
% 
% subplot(2,2,3.5)
% plot(Mo(1,2:end),Mo(2,2:end),'.m');hold on;
% plot(Sep(1,:),Sep(2,:),'.k','markersize',25);axis([0 4 -1.5 2.0])
% for k = 1:length(Sep)
%    text(Sep(1,k)+0.1,Sep(2,k)+0.1,[num2str(k)],'fontsize',14) ;
% end
% set(gca,'FontName','Times','FontSize',18);
% set(gca,'TickLabelInterpreter','latex');
% xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",16);
% ylabel('$Z\ (\rm m)$','Interpreter','latex',"FontSize",16);
% box on;grid on;