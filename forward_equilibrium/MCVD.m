function [PSI] = MCVD(rr1,zz1,CRD,CZD,CJ)
% MCVD : 计算线圈的磁通量 update:更具输入Coils返回线圈生成的矩阵
% 1. 当R中包含0时，将psi(0,:) = 0;
% 2. 当R中没有0时，使用椭圆积分计算psi
c0=2e-7;
[RR,ZZ]=ndgrid(rr1,zz1);
PSI = zeros(size(RR));

NR = size(CRD,1);NZ = size(CRD,2);
I0 = CJ;

if ( rr1(1) == 0 )
    tmpRR = RR(2:end,:);tmpZZ = ZZ(2:end,:);
    Itmp=I0/(NR*NZ);
    
    for jr=1:NR
        for jz=1:NZ
            tmpCR = CRD(jr,jz);
            tmpCZ = CZD(jr,jz);
            k2 = 4*tmpCR*tmpRR./((tmpRR + tmpCR).^2 + (tmpZZ - tmpCZ).^2);
            [MK,ME] = ellipke(k2);
            Psi_tmp=-c0*Itmp.*tmpRR.*sqrt(tmpCR./tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
            PSI(2:end,:)=PSI(2:end,:) + Psi_tmp;
        end
    end
    
    PSI(1,:)=0;
    
else
    tmpRR = RR;tmpZZ = ZZ;
    Itmp=I0/(NR*NZ);
    for jr=1:NR
        for jz=1:NZ
            tmpCR = CRD(jr,jz);
            tmpCZ = CZD(jr,jz);
            k2 = 4*tmpCR*tmpRR./((tmpRR + tmpCR).^2 + (tmpZZ - tmpCZ).^2);
            [MK,ME] = ellipke(k2);
            Psi_tmp = -c0*Itmp.*tmpRR.*sqrt(tmpCR./tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
            PSI = PSI + Psi_tmp;
        end
    end
    
end
% 在线圈位置位于网格上时，psi表达式无定义，使用临近点插值计算无定义点
[row,col] = find(PSI == -inf);
for j=1:length(row)
    PSI(row(j),col(j))=(PSI(row(j)+1,col(j)+1)+PSI(row(j)+1,col(j)+0)+...
        PSI(row(j)+1,col(j)-1)+PSI(row(j)-1,col(j)+1)+...
        PSI(row(j)-1,col(j)+0)+PSI(row(j)-1,col(j)-1)+...
        PSI(row(j),col(j)+1)+PSI(row(j),col(j)-1))/8;
end
