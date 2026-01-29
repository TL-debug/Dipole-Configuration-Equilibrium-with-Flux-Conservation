function [NCD,MuT] = SDP(RR,ZZ,VCR,VCZ)
% 计算线圈对计算域的互感
mu0 = 4*pi*10^-7;
[CRD,~] = ndgrid(VCR,VCZ);
NCD = size(CRD)*flipud(size(CRD)')/2;%上下倒置
JRD = repmat(VCR,1,length(VCZ));
JZD = reshape(repmat(VCZ,length(VCR),1),[],1)';
MuT = zeros(length(RR(:,1)),length(ZZ(1,:)),NCD);
for j = 1:NCD
    k2 = 4*JRD(j).*RR(2:end-1,2:end-1)./...
        ( (RR(2:end-1,2:end-1) + JRD(j)).^2 + (ZZ(2:end-1,2:end-1) - JZD(j)).^2 );
    [MK,ME] = ellipke(k2);
    MuT(2:end-1,2:end-1,j) = mu0*sqrt(RR(2:end-1,2:end-1).*JRD(j))./...
                             sqrt(k2).*((2 - k2).*MK - 2*ME);
end
MuT(find(abs( MuT ) > 10 )) = 0;
end