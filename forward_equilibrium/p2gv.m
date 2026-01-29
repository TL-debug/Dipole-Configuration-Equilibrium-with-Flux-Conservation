function [RHSC] = p2gv(Coils,GR,GZ)
%	update:2024/10/19 test ok!
%   P2GV 按照输入线圈的顺序，将电流值由线性插值贡献给相邻周围的四个网格
%   输入：Cell线圈的位置和对应线圈电流信息 ; GR,GZ 网格的向量
%       1.Cell矩阵第一列是线圈的径向位置，第二列是线圈的纵向位置，...
%         第三列为线圈的电流，每一行表示线圈编号。
%       2.支持矩阵输入
%   输出 RHSC 

% 循环对RHSC赋值，且不计算Coils中计算域外的线圈电流对RHSC的贡献
rr = GR;dr = rr(2) - rr(1);
zz = GZ;dz = zz(2) - zz(1);
RHSC = zeros(length(rr),length(zz));
for k = 1:size(Coils,1)
    RHSCt = 0*RHSC;
    % 读取线圈k的位置和电流信息
    CR = reshape(cell2mat(Coils(k,1)),[],1);
    CZ = reshape(cell2mat(Coils(k,2)),[],1);
    CI = cell2mat(Coils(k,3))/length(CR);
    
    % 判断是否在计算域内(在计算域内则将电流贡献分配到临近的网格)
    if ( max(zz) > max(CZ) )
        % R 方向指标
        xrow = floor( (CR - rr(1))./dr ) + 1;  % 点位置
        xcol = 1:length(CR);
        dxf  = CR - ( ( xrow - 1 )*dr + rr(1) );
        dxb  = dr - dxf;
        
        % Z 方向指标
        yrow = floor( (CZ - zz(1))./dz ) + 1;
        ycol = 1:length(CZ);
        dyf = CZ - ( ( yrow - 1 )*dz + zz(1) );
        dyb = dz - dyf;
        
        MX = sparse(xrow+1,xcol',dxf,length(rr),length(CR))+...
            sparse(xrow,xcol',dxb,length(rr),length(CR));
        MY = sparse(yrow+1,ycol',dyf,length(zz),length(CZ))+...
            sparse(yrow,ycol',dyb,length(zz),length(CZ));
        MY = MY';
        I = CI*speye(length(CZ));
        RHSCt = MX*I*MY/(dr*dz);
        
    end
    RHSC = RHSC + RHSCt;
end

end

