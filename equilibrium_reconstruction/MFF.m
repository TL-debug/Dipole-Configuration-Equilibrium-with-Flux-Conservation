function [PSI] = MFF(rr1,zz1,Coils)
% MFF magnetic flux of current filaments  计算线圈的磁通量
%   输入Coils R，Z
%   功能 1. 当R中包含0时，将psi(0,:) = 0;
%        2. 当R中没有0时，使用椭圆积分计算psi
c0=2e-7;
[RR,ZZ]=ndgrid(rr1,zz1);
PSI = zeros(length(rr1),length(zz1),size(Coils,1));

for k1 = 1:size(Coils,1)
    CRD = cell2mat(Coils(k1,1));
    CZD = cell2mat(Coils(k1,2));
    NR = size(CRD,1);
    NZ = size(CRD,2);
    I0 = cell2mat(Coils(k1,3));
    tPSI = 0*RR;
    if ( rr1(1) == 0 )
        tmpRR = RR(2:end,:);tmpZZ = ZZ(2:end,:);
        Itmp=I0/(NR*NZ);
        
        for jr=1:NR
            for jz=1:NZ
                tmpCR = CRD(jr,jz);
                tmpCZ = CZD(jr,jz);
                k2 = 4*tmpCR*tmpRR./((tmpRR + tmpCR).^2 + (tmpZZ - tmpCZ).^2);
                [MK,ME] = ellipke(k2);
                Psi_tmp= -c0*Itmp.*tmpRR.*sqrt(tmpCR./tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
                tPSI(2:end,:) = tPSI(2:end,:) + Psi_tmp;
            end
        end

        tPSI(1,:)=0;

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
                tPSI = tPSI + Psi_tmp;
            end
        end

    end
    % 在线圈位置位于网格上时，psi表达式无定义，使用临近点插值计算无定义点
    [row,col] = find(abs(tPSI) == inf);
    for j=1:length(row)
        tPSI(row(j),col(j))=(tPSI(row(j)+1,col(j)+1)+tPSI(row(j)+1,col(j)+0)+...
            tPSI(row(j)+1,col(j)-1)+tPSI(row(j)-1,col(j)+1)+...
            tPSI(row(j)-1,col(j)+0)+tPSI(row(j)-1,col(j)-1)+...
            tPSI(row(j),col(j)+1)+tPSI(row(j),col(j)-1))/8;
    end

    PSI(:,:,k1) = tPSI; 
end

end