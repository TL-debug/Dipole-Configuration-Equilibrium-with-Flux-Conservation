function [NCD,MuT] = GMI(RR,ZZ,CRD,CZD)
% GMI 线圈与网格之间的互感   update：2024/10/22
%   输入网格与线圈的位置信息；输出每个线圈丝与计算域网格的互感矩阵
mu0 = 4*pi*10^-7;
VRD = reshape(CRD,[],1);
VZD = reshape(CZD,[],1);
NCD = length(VRD);
MuT = zeros(length(RR(:,1)),length(ZZ(1,:)),NCD);
tRR = RR(2:end-1,2:end-1);
tZZ = ZZ(2:end-1,2:end-1);
for k1 = 1:NCD
    % k2不计算R=0的点，
    k2 = 4*VRD(k1)*tRR./( (tRR + VRD(k1)).^2 + (tZZ - VZD(k1)).^2 );
    [MK,ME] = ellipke(k2);
    tmp_MuT = mu0*sqrt(tRR.*VRD(k1)./k2).*( (2 - k2).*MK - 2*ME );
    if ( max ( tmp_MuT == Inf ) == 1 )
        % - 插值计算无穷点
        ind_Minf = find( tmp_MuT == inf );
        mut_int  = interp2(tRR',tZZ',tmp_MuT',tRR(ind_Minf),tZZ(ind_Minf));
        tmp_MuT ( ind_Minf ) = mut_int;
    end
    MuT(2:end-1,2:end-1,k1) = tmp_MuT;
end

