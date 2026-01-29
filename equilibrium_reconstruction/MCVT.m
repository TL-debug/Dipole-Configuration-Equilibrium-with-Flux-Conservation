function [PSI] = MCVT(RR,ZZ,Coils)
% MCVT 通用的线圈通量生成子程序
% 输入：  网格R和网格Z；线圈位置和电流的cell向量第一列R；第二列Z；第三列为电流
% Note:  (不管匝数大小都是用统一的cell保存）
% 输入：磁通量
c0  = 2e-7;
PSI = zeros(size(RR));
dr  = RR(2,1);dz = ZZ(1,2)-ZZ(1,1);

for j1 = 1:size(Coils,1)
    tmpPSI = zeros(size(RR));
    CRD  = cell2mat(Coils(j1,1));
    CZD  = cell2mat(Coils(j1,2));
    NR = size(CRD,1);
    NZ = size(CRD,2);
    Itmp   =  cell2mat(Coils(j1,3))/(NR*NZ);

    if ( RR(1) == 0 )
        tmpRR = RR(2:end,:);tmpZZ = ZZ(2:end,:);

        for jr=1:NR
            for jz=1:NZ
                tmpCR = CRD(jr,jz);
                tmpCZ = CZD(jr,jz);
                k2 = 4*tmpCR*tmpRR./((tmpRR + tmpCR).^2 + (tmpZZ - tmpCZ).^2);
                [MK,ME] = ellipke(k2);
                Psi_tmp=-c0*Itmp.*tmpRR.*sqrt(tmpCR./tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
                
                 % 找到无定义点，并使用dr和dz进行赋值
                [row,col] = find(Psi_tmp == -inf | Psi_tmp == inf );

                for jnp = 1:length(row)% 网格位置
                    tmpRrlud = [tmpRR(row(jnp),col(jnp)) + dr/16,tmpRR(row(jnp),col(jnp))-dr/16,...
                                tmpRR(row(jnp),col(jnp)),tmpRR(row(jnp),col(jnp))];
                    tmpZrlud = [tmpZZ(row(jnp),col(jnp)),tmpZZ(row(jnp),col(jnp)),...
                                tmpZZ(row(jnp),col(jnp)) + dr/16,tmpZZ(row(jnp),col(jnp))-dr/16];
                    k2rlud   = 4*tmpRrlud*tmpCR./( (tmpRrlud + tmpCR).^2 + (tmpZrlud - tmpCZ).^2 );
                    [Krlud,Erlud] =  ellipke(k2rlud);
                    tmpGrlud = -c0*Itmp.*sqrt(tmpRrlud*tmpCR)./sqrt(k2rlud).*((2-k2rlud).*Krlud - 2*Erlud);

                    Psi_tmp(row(jnp),col(jnp)) = sum(tmpGrlud)/length(tmpGrlud);
                    
                end
                   
                PSI(2:end,:)=PSI(2:end,:) + Psi_tmp;

            end
        end

        PSI(1,:)=0;

    else
        tmpRR = RR;tmpZZ = ZZ;

        for jr=1:NR
            for jz=1:NZ
                tmpCR = CRD(jr,jz);
                tmpCZ = CZD(jr,jz);
                k2 = 4*tmpCR*tmpRR./((tmpRR + tmpCR).^2 + (tmpZZ - tmpCZ).^2);
                [MK,ME] = ellipke(k2);
                Psi_tmp = -c0*Itmp.*sqrt(tmpCR.*tmpRR)./sqrt(k2).*((2-k2).*MK-2*ME);
                
                
                % 找到无定义点，并使用dr和dz进行赋值
                [row,col] = find(Psi_tmp == -inf | Psi_tmp == inf );

                for jnp = 1:length(row)% 网格位置
                    tmpRrlud = [RR(row(jnp),col(jnp))+dr/16,RR(row(jnp),col(jnp))-dr/16,...
                                RR(row(jnp),col(jnp)),RR(row(jnp),col(jnp))];
                    tmpZrlud = [ZZ(row(jnp),col(jnp)),ZZ(row(jnp),col(jnp)),...
                                ZZ(row(jnp),col(jnp)) + dr/16,ZZ(row(jnp),col(jnp))-dr/16];
                    k2rlud   = 4*tmpRrlud*tmpCR./( (tmpRrlud + tmpCR).^2 + (tmpZrlud - tmpCZ).^2 );
                    [Krlud,Erlud] =  ellipke(k2rlud);
                    tmpGrlud = -c0*Itmp.*sqrt(tmpRrlud*tmpCR)./sqrt(k2rlud).*((2-k2rlud).*Krlud - 2*Erlud);


                    Psi_tmp(row(jnp),col(jnp)) = sum(tmpGrlud)/length(tmpGrlud);
                    
                end
                
                PSI = PSI + Psi_tmp;

            end
        end

    end
    PSI = PSI + tmpPSI;
end

% % 在线圈位置位于网格上时，psi表达式无定义，使用临近点插值计算无定义点
% [row,col] = find(PSI == -inf | PSI == inf ); % 四点drdz进行计算无定于点
% for j=1:length(row)
%     PSI(row(j),col(j))=(PSI(row(j)+1,col(j)+1)+PSI(row(j)+1,col(j)+0)+...
%         PSI(row(j)+1,col(j)-1)+PSI(row(j)-1,col(j)+1)+...
%         PSI(row(j)-1,col(j)+0)+PSI(row(j)-1,col(j)-1)+...
%         PSI(row(j),col(j)+1)+PSI(row(j),col(j)-1))/8;
% end

end

