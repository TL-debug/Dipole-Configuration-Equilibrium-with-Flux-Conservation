function [MD] = FSI(CRD,CZD)
% Filament Method compute Self-Introduce
mu0 = 4*pi*10^-7;
VCR = CRD(:,1);
VCZ = CZD(1,:);
MutD = 0*CRD;
for i = 1:length(VCR)
    for j = 1:length(VCZ)
        % - ignore 
        k2 = 4*VCR(i).*CRD./((VCR(i) + CRD).^2 + (VCZ(j) - CZD).^2);
        [MK,ME] = ellipke(k2);
        tMutD = mu0*sqrt(CRD.*VCR(i)./k2).*((2 - k2).*MK - 2*ME);
        tMutD(i,j) = 0;
        MutD(i,j) = sum(sum(tMutD));
    end
end
MD = sum(sum(MutD));
end

