function [Psi] = MGSOR2(Psi_S1,RES,MRI,NIT,ome)
% - Update:2024/7/1
% MRI 第一列弛豫矩阵 第二列插值矩阵 第三和第四列表示网格
% ---- V Cycle 测试ok 
% EV = {};%用于记录迭代后的结果作为初值
% for k = 1:2*size(MRI,1) - 1 % 每次迭代传递残差
%     if ( k <= size(MRI,1) )   % 网格向下（粗）% 保留迭代结果
%         MR  = cell2mat(MRI(k,1));
%         RHS = MR*RES(:);RHS = reshape(RHS,sqrt(size(MR,1)),[]);
%         E0  = zeros(size(RHS));
%         MR  = cell2mat(MRI(k,3));  
%         MZ  = cell2mat(MRI(k,4)); 
%         [E1,RES] = GSSOR2M(E0,RHS,MR,MZ,NIT,ome);
%         EV = [EV;{E1}];
%         if ( k == size(MRI,1) )
%             E0D = E1;
%         end
%     elseif ( k > size(MRI,1) )% 网格向上（细）%  迭代结果向上传递
%         MI  = cell2mat(MRI(2*size(MRI,1) - k + 1,2));
%         RHS = MI*RES(:); RHS = reshape(RHS,sqrt(size(MI,1)),[]);
%         
%         E0  = reshape(MI*E0D(:),sqrt(size(MI,1)),[]) +...
%               cell2mat(EV(2*size(MRI,1) - k));
%           
%         MR  = cell2mat(MRI(2*size(MRI,1) - k,3));
%         MZ  = cell2mat(MRI(2*size(MRI,1) - k,4));
%         [E1,RES] = GSSOR2M(E0,RHS,MR,MZ,NIT,ome);
%         E0D = E1;
%     end  
% end
%     Psi = Psi_S1 + reshape(cell2mat(MRI(1,2))*E0D(:),sqrt(size(cell2mat(MRI(1,2)),1)),[]);
% end

% ---- Full V Cycle 速度最慢
% - 计算每层的RHS(矩阵)项保存
RHS = [];
for k = 1:size(MRI,1)
    tmpMR = cell2mat(MRI(k,1));
    if ( k == 1 )
        tmpRHS = reshape(tmpMR*RES(:),sqrt(size(tmpMR,1)),[]);
    else
        tRHS   = cell2mat(RHS(k-1));
        tmpRHS = reshape(tmpMR*tRHS(:),sqrt(size(tmpMR,1)),[]);
    end
    RHS = [RHS;{tmpRHS}];
end

indk = flip(1:length(RHS));
for k = 1:length(indk)
    % - 进行最底层迭代，传递迭代值
    if ( k == 1 )
        rhs  = cell2mat(RHS(indk(k)));
        mr   = cell2mat(MRI(indk(k),3));
        mz   = cell2mat(MRI(indk(k),4));
        [E1,~] = GSSOR2M(zeros(size(rhs)),rhs,mr,mz,NIT,ome);
        % - 进行向上的V循环
    else
        [E1] = MGFV(MRI(indk(k):end,:),RHS(indk(k)),E0,NIT,ome);
    end
    % - 插值到上一层作为初值
    E0 = reshape(cell2mat(MRI(indk(k),2))*E1(:),sqrt(size(cell2mat(MRI(indk(k),2)),1)),[]);
end

Psi = Psi_S1 + E0;

end

% ---- W Cycle waiting add .....


% ---- N Cycle Waiting add ......




% 	gam=zeros(size(RES));
%     % - F cycle 将RES4 返回上一步进行迭代 进行类似的迭代
%     % - 从最粗的网格开始
%     % - 第一层
%     [~,~,MRB,MIB] = RB2D(RES); % MRB 由粗到细
% 	  sz=size(MRB);
%     RHS1 = MRB*RES(:);RHS1 = reshape(RHS1,sqrt(sz(1)),[]);%残差粗化2
%     R1= MRB*R(:);R1 = reshape(R1,sqrt(sz(1)),[]);
%     Z1= MRB*Z(:);Z1 = reshape(Z1,sqrt(sz(1)),[]);
% 
%     
%     % - 第二层
% 	 [~,~,MRB1,MIB1] = RB2D(RHS1);% MRB 由粗到细
%     sz=size(MRB1);
%     RHS2 = MRB1*RHS1(:);RHS2 = reshape(RHS2,sqrt(sz(1)),[]);%残差粗化3
%     R2 = MRB1*R1(:);R2 = reshape(R2,sqrt(sz(1)),[]);
%     Z2 = MRB1*Z1(:);Z2 = reshape(Z2,sqrt(sz(1)),[]);
% 
%     
%     % - 第三层
% 	[~,~,MRB2,MIB2] = RB2D(RHS2);
%     sz = size(MRB2);
%     RHS3 = MRB2*RHS2(:);RHS3 = reshape(RHS3,sqrt(sz(1)),[]);%残差粗化4
%     R3 = MRB2*R2(:);R3 = reshape(R3,sqrt(sz(1)),[]);
%     Z3 = MRB2*Z2(:);Z3 = reshape(Z3,sqrt(sz(1)),[]);
% 
% 
% FMG Cycle
    
    
    % ---------- First Layer ---------- % % 该多网格方法没有迭代返回的参数
% 	[~,~,MRB,MIB] = RB2D(RES);
%     sz=size(MRB);
%     RHS1=MRB*RES(:);RHS1=reshape(RHS1,sqrt(sz(1)),[]);
%     R1=MRB*R(:);R1=reshape(R1,sqrt(sz(1)),[]);
%     Z1=MRB*Z(:);Z1=reshape(Z1,sqrt(sz(1)),[]);
%     E01=zeros(size(RHS1));
%     [E10,RES1] = GSSOR2M(E01,RHS1,R1,Z1,NIT,ome);
% %     gamt = MIB*E10(:);
% %     gam = gam + reshape(gamt,sqrt(length(gamt)),[]);
%     
%     [~,~,MRB1,MIB1] = RB2D(RHS1);
%     sz=size(MRB1);
%     RHS2=MRB1*RES1(:);RHS2=reshape(RHS2,sqrt(sz(1)),[]);
%     R2=MRB1*R1(:);R2=reshape(R2,sqrt(sz(1)),[]);
%     Z2=MRB1*Z1(:);Z2=reshape(Z2,sqrt(sz(1)),[]);
%     E02=zeros(size(RHS2));
%     [E20,RES2] = GSSOR2M(E02,RHS2,R2,Z2,NIT,ome);
% %     gamt=MIB*MIB1*E20(:);
% %     gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
%     [~,~,MRB2,MIB2] = RB2D(RHS2);
%     sz=size(MRB2);
%     RHS3=MRB2*RES2(:);RHS3=reshape(RHS3,sqrt(sz(1)),[]);
%     R3=MRB2*R2(:);R3=reshape(R3,sqrt(sz(1)),[]);
%     Z3=MRB2*Z2(:);Z3=reshape(Z3,sqrt(sz(1)),[]);
% 	E03=zeros(size(RHS3));
%     [E30,RES3] = GSSOR2M(E03,RHS3,R3,Z3,NIT,ome);
% % 	gamt=MIB*MIB1*MIB2*E30(:);
% %     gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
%     [~,~,MRB3,MIB3] = RB2D(RHS3);
%     sz=size(MRB3);
%     RHS4=MRB3*RES3(:);RHS4=reshape(RHS4,sqrt(sz(1)),[]);
%     R4=MRB3*R3(:);R4=reshape(R4,sqrt(sz(1)),[]);
%     Z4=MRB3*Z3(:);Z4=reshape(Z4,sqrt(sz(1)),[]);
% 	E04=zeros(size(RHS4));
% 	[E40,RES4] = GSSOR2M(E04,RHS4,R4,Z4,NIT,ome);
% %     gamt=MIB*MIB1*MIB2*MIB3*E40(:);
% %     gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
% 
%     [~,~,MRB4,MIB4] = RB2D(RHS4);
% 	sz=size(MRB4);
%     RHS5=MRB4*RES4(:);RHS5=reshape(RHS5,sqrt(sz(1)),[]);
%     R5=MRB4*R4(:);R5=reshape(R5,sqrt(sz(1)),[]);
%     Z5=MRB4*Z4(:);Z5=reshape(Z5,sqrt(sz(1)),[]);
% 	E05=zeros(size(RHS5));
%     [E50,RES5] = GSSOR2M(E05,RHS5,R5,Z5,NIT,ome);
% % 	gamt=MIB*MIB1*MIB2*MIB3*MIB4*E50(:);
% %     gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
%     [~,~,MRB5,MIB5] = RB2D(RHS5);
% 	sz=size(MRB5);
%     RHS6=MRB5*RES5(:);RHS6=reshape(RHS6,sqrt(sz(1)),[]);
%     R6=MRB5*R5(:);R6=reshape(R6,sqrt(sz(1)),[]);
%     Z6=MRB5*Z5(:);Z6=reshape(Z6,sqrt(sz(1)),[]);
% 	E06=zeros(size(RHS6));
% 	[E60,RES6] = GSSOR2M(E06,RHS6,R6,Z6,NIT,ome);
% % 	gamt=MIB*MIB1*MIB2*MIB3*MIB4*MIB5*E60(:);
% %   gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
% 
%     gamt = MIB5*E60(:);
% 	[tE50,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + E50,RES5,R5,Z5,NIT,ome);
%     gamt = MIB4*tE50(:);
%     [tE40,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + E40,RES4,R4,Z4,NIT,ome);
%     gamt = MIB3*tE40(:);
%     [tE30,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + E30,RES3,R3,Z3,NIT,ome);
%     gamt = MIB2*tE30(:);
%     [tE20,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + E20,RES2,R2,Z2,NIT,ome);
%     gamt = MIB1*tE20(:);
%     [tE10,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + E10,RES1,R1,Z1,NIT,ome);
%     gamt = MIB*tE10(:);
%     [tE00,~] = GSSOR2M(reshape(gamt,sqrt(length(gamt)),[]) + Psi_S1,RES,R,Z,NIT,ome);
    

    
%     [~,~,MRB6,MIB6] = RB2D(RHS6);
%     sz=size(MRB6);
%     RHS7=MRB6*RES6(:);RHS7=reshape(RHS7,sqrt(sz(1)),[]);
%     R7=MRB6*R6(:);R7=reshape(R7,sqrt(sz(1)),[]);
%     Z7=MRB6*Z6(:);Z7=reshape(Z7,sqrt(sz(1)),[]);
%     E00=zeros(size(RHS7));
%     [E70,~] = GSSOR2M(E00,RHS7,R7,Z7,NIT,ome);
%     gamt=MIB*MIB1*MIB2*MIB3*MIB4*MIB5*MIB6*E70(:);
%     gam=gam+reshape(gamt,sqrt(length(gamt)),[]);
%     
%     E0 = Psi_S1 + gam;  % 这里残差的边界值为0