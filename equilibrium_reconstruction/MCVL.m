function [PSI] = MCVL(rr1,zz1,CRL,CZL,CJ)
c0=2e-7;
[RR,ZZ]=ndgrid(rr1,zz1);
PSI=zeros(size(RR));Psi_tmp = PSI;

tmp1=size(CRL);
NR=tmp1(1);NZ=tmp1(2);
I0=CJ;
Itmp=I0/(NR*NZ);
tmpRR=RR(2:end,:);tmpZZ=ZZ(2:end,:);
    for jr=1:NR
        for jz=1:NZ
            tmpCR=CRL(jr,jz);
            tmpCZ=CZL(jr,jz);
            k2=4*tmpCR*tmpRR./((tmpRR+tmpCR).^2+(tmpZZ-tmpCZ).^2);
            [MK,ME] = ellipke(k2);
            Psi_tmp=-c0*Itmp.*tmpRR.*sqrt(tmpCR./tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
            PSI(2:end,:)=PSI(2:end,:)+Psi_tmp;
        end
    end
PSI(1,:)=0;
[row,col] = find(PSI == -inf);
    for j=1:length(row)
        PSI(row(j),col(j))=(PSI(row(j)+1,col(j)+1)+PSI(row(j)+1,col(j)+0)+...
                           PSI(row(j)+1,col(j)-1)+PSI(row(j)-1,col(j)+1)+...
                           PSI(row(j)-1,col(j)+0)+PSI(row(j)-1,col(j)-1)+...
                           PSI(row(j),col(j)+1)+PSI(row(j),col(j)-1))/8;
    end
