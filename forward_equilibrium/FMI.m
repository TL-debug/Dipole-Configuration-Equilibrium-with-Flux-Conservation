function [mDL] = FMI(CRD,CZD,CRL,CZL) % dimensions of the D and L coils
% FMI: Filament Method Multi-Introduce
mu0 = 4*pi*10^-7;
MDL = 0*CRD;
trr = CRD(:,1);
tzz = CZD(1,:);
for i = 1:length(trr)
    for j = 1:length(tzz)
        k2 = 4*trr(i).*CRL./((trr(i) + CRL).^2 + (tzz(j) - CZL).^2);
        [MK,ME] = ellipke(k2);
        tMDL = mu0*sqrt(trr(i).*CRL./k2).*((2 - k2).*MK - 2*ME);
        MDL(i,j) = sum(sum(tMDL));
    end
end
mDL = sum(sum(MDL));
end

