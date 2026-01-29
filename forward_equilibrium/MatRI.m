function [MIR] = MatRI(R,Z,nlay)
TR = R;
TZ = Z;
MIR = {};
for k = 1:nlay
    [~,~,tMRB,tMIB] = RB2D(TR);
    sz = size(tMRB,1);
    tR = tMRB*TR(:);tR = reshape(tR,sqrt(sz),[]);
    tZ = tMRB*TZ(:);tZ = reshape(tZ,sqrt(sz),[]);
    
    RIRZ = [{tMRB},{tMIB},{tR},{tZ}];
    MIR = [MIR;RIRZ];
    
    TR = tR;
    TZ = tZ;
end
end

% % 测试代码 %
% % initial Sys
% % nI=9;nlay=7;nbl=2^(nI-nlay);
% % NG=2^nlay*nbl-(2^nlay-1); % 包括边界上网格
% nI = 7;NG = 2^nI + 1;nlay = nI - 3;%多重网格层数
% NR=NG;NZ=NG;c0=2e-7;r1=0;r2=4;z1=-1.5;z2=2.50;mu0=4*pi*1e-7;
% rr=linspace(r1,r2,NR);zz=linspace(z1,z2,NZ);dr=rr(2)-rr(1);dz=zz(2)-zz(1);
% %------生成空间中的磁场和通量子程序
% [R,Z]=ndgrid(rr,zz);

% % ----- 返回矩阵形式
% tic
% for j = 1:10000
%     
%     TR = R;
%     TZ = Z;
%     MIR = {};
%     for k = 1:nlay
%         [~,~,tMRB,tMIB] = RB2D(TR);
%         sz = size(tMRB,1);
%         tR = tMRB*TR(:);tR = reshape(tR,sqrt(sz),[]);
%         tZ = tMRB*TZ(:);tZ = reshape(tZ,sqrt(sz),[]);
%         
%         RIRZ = [{tMRB},{tMIB},{tR},{tZ}];
%         MIR = [MIR;RIRZ];
%         
%         TR = tR;
%         TZ = tZ;
%     end
%     j = j + 1;
% end
% toc
%历时 10.823579 秒。