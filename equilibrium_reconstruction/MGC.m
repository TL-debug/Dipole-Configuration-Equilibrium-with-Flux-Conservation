function [NCD,MuT] = MGC(RR,ZZ,CRD,CZD)
% MGC 线圈和网格之间的互感
% 输入：网格与网格的相对位置
% 输出：MuT 行表示线圈；列表示网格(respahe(tmpJ1,[],1))
%      ncd D输入线圈的丝电流个数
mu0 = 4*pi*10^-7;
NCD = length(reshape(CRD,[],1));
MuT = zeros(NCD,length(reshape(RR,[],1)));


VRD = reshape(CRD,[],1);
VZD = reshape(CZD,[],1);

% Note:
% 1.当R=0时不计算R=0网格与线圈之间的互感；R<0时则计算所有网格与线圈之间的互感
% 2.删除网格与线圈重复位置的通量，并插值重写计算该点的互感值

% R包含0的网格与线圈之间的互感;先计算再补全然后再reshape
tR = RR(2:end-1,2:end-1);
tZ = ZZ(2:end-1,2:end-1);
tM = 0*RR;
for k1 = 1:NCD
    k2 = 4*tR*VRD(k1)./( (tR + VRD(k1)).^2 + (tZ - VZD(k1)).^2 );
    [MK,ME] = ellipke(k2);
    tmp_MuT = mu0*sqrt(tR.*VRD(k1)./k2).*( (2 - k2).*MK - 2*ME );
    % 删除网格与线圈重复位置的通量，并插值重写计算该点的互感值
    if ( max ( tmp_MuT == Inf ) == 1 )
        % - 插值计算无穷点
        ind_Minf = find( tmp_MuT == inf );
        mut_int  = interp2(tR',tZ',tmp_MuT',tR(ind_Minf),tZ(ind_Minf));
        tmp_MuT ( ind_Minf ) = mut_int;
    end
    tM(2:end-1,2:end-1) = tmp_MuT;
    MuT(k1,:) = reshape(tM,[],1)';
end

% elseif( RR(1,1) ~=0 )
%     % R不包含0的网格与线圈之间的互感
%     tR = RR;
%     tZ = ZZ;
%     for k1 = 1:NCD
%         k2 = 4*tR*VRD(k1)./( (tR + VRD(k1)).^2 + (tZ - VZD(k1)).^2 );
%         [MK,ME] = ellipke(k2);
%         tmp_MuT = mu0*sqrt(tR.*VRD(k1)./k2).*( (2 - k2).*MK - 2*ME );
%         % 删除网格与线圈重复位置的通量，并插值重写计算该点的互感值
%         if ( max ( tmp_MuT == Inf ) == 1 )
%             % - 插值计算无穷点
%             ind_Minf = find( tmp_MuT == inf );
%             mut_int  = interp2(tR',tZ',tmp_MuT',tR(ind_Minf),tZ(ind_Minf));
%             tmp_MuT ( ind_Minf ) = mut_int;
%         end
%         MuT(k1,:) = reshape(tmp_MuT,[],1);
%     end
% end
end

% % ------------- 参考公式
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