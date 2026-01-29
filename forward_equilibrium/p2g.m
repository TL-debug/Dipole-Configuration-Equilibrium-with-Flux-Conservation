function [RHS] = p2g(LDR,LDZ,LLR,LLZ,IDL,GR,GZ)
% p2g points to grids 点的数值由线性插值贡献给相邻周围的四个网格
% ---- 更新p2g update：2024/10/2 
% 一个更为稳定的点数值插值到网格程序,避免原程序中由于R未从0取导致的RHSC分布存在复数的问题
% 当计算网格中包括0时，p2g的结果不正确的
    xx = GR;dx = xx(2) - xx(1);
    yy = GZ;dy = yy(2) - yy(1);
    % ---------- 线圈位置和电流初始化 ---------- %
    LDX = reshape(LDR,[],1);
    LDY = reshape(LDZ,[],1);
    CDI = IDL(2)/length(LDX); % 偶极场线圈
    LLX = reshape(LLR,[],1);
    LLY = reshape(LLZ,[],1);
    CLI = IDL(1)/length(LLX); % 悬浮线圈


    % R 方向指标
    xrow = floor( (LDX - xx(1))./dx ) + 1;  % 点位置
    xcol = 1:length(LDX);                   % 
    
    dxf = LDX - ( ( xrow - 1 )*dx + xx(1) );% 对前一个网格贡献
%     dxf = LDX - (floor((LDX-xx(1))./dx))*dx;% 对前一个网格贡献
    dxb = dx - dxf;                         % 对后一个网格贡献
%     plot(dxf,'.b');hold on;plot(dxb,'.r');
    % Z 方向指标
    yrow = floor( (LDY - yy(1))./dy )+1;
    ycol = 1:length(LDY);
    
    dyf = LDY - ( ( yrow - 1 )*dy + yy(1) );
%     dyf=LDY-floor(LDY./dy)*dy;
    dyb = dy - dyf;
    
    MX=sparse(xrow+1,xcol',dxf,length(xx),length(LDX))+...
       sparse(xrow,xcol',dxb,length(xx),length(LDX));
    MY=sparse(yrow+1,ycol',dyf,length(yy),length(LDY))+...
       sparse(yrow,ycol',dyb,length(yy),length(LDY));
    MY=MY';
    I=CDI*speye(length(LDY));
    RHS = MX*I*MY/(dx*dy);
    
    if ( max(yy) >=  max(LLY) )  %求解区域包含偶极场和悬浮线圈 
        xrow = floor( (LLX - xx(1))./dx ) + 1;
        xcol = 1:length(LLX);
        dxf  = LLX - ( ( xrow -1 )*dx + xx(1) );
        dxb  = dx - dxf;
        % Z 方向指标
        yrow = floor( (LLY - yy(1))./dy ) + 1;
        ycol = 1:length(LLY);
        dyf  = LLY - ( ( yrow - 1 )*dy + yy(1) );
        dyb  = dy-dyf;
        
        MX=sparse(xrow+1,xcol',dxf,length(xx),length(LLX))+...
           sparse(xrow,xcol',dxb,length(xx),length(LLX));
        MY=sparse(yrow+1,ycol',dyf,length(yy),length(LLY))+...
           sparse(yrow,ycol',dyb,length(yy),length(LLY));
        MY=MY';
        I=CLI*speye(length(LLY));
        RHS1 = MX*I*MY/(dx*dy);
        RHS = RHS + RHS1;
    end
    RHS = full(RHS);