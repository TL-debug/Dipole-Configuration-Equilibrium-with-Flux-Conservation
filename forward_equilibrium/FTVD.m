function [Vf,Dvdp,Dpdp] = FTVD(psiout,tmpP1,BBR,BBZ,R,Z,indent,indSepz)
% FTVD 通量管体积以及对R和Z的偏导
% 输入： 通量 磁场 R，Z，参考点； 输出：通关体积，通量体积对通量导数，压强对通量的导数（链式法则）

Dvdp = 0*indent;
Vf   = 0*indent;
Dpdp = 0*indent;

dr   = R(2,1)-R(1,1);
dz   = Z(1,2)-Z(1,1);
for jr = 1:length(indent)
% - 一个对应的参考点需要计算5个通量体积
    % 确定包括参考点的9点位置R方向5点；Z方向4
    rind = indent(jr);zind = indSepz;
%     tmpvr = [R(rind-2,zind) R(rind-1,zind) R(rind,zind) R(rind+1,zind) R(rind+2,zind)... % R 方向
%              R(rind,zind-2) R(rind,zind-1) R(rind,zind+1) R(rind,zind+2)];               % Z 方向
%     tmpvz = [Z(rind-2,zind) Z(rind-1,zind) Z(rind,zind) Z(rind+1,zind) Z(rind+2,zind)... % R 方向
%              Z(rind,zind-2) Z(rind,zind-1) Z(rind,zind+1) Z(rind,zind+2)];               % Z 方向
    tmpvp = [psiout(rind-2,zind) psiout(rind-1,zind) psiout(rind,zind)...    % 参考点的通量
             psiout(rind+1,zind) psiout(rind+2,zind) psiout(rind,zind-2)...
             psiout(rind,zind-1) psiout(rind,zind+1) psiout(rind,zind+2)];

    tmpvr = [tmpP1(rind-2,zind) tmpP1(rind-1,zind) tmpP1(rind,zind)...       % 参考点的压强
             tmpP1(rind+1,zind) tmpP1(rind+2,zind) tmpP1(rind,zind-2)...
             tmpP1(rind,zind-1) tmpP1(rind,zind+1) tmpP1(rind,zind+2)];
    

% 通量对R和Z的导数
    dpdr = (tmpvp(1) - 8*tmpvp(2) + 8*tmpvp(4) - tmpvp(5))/(12*dr);
    dpdz = (tmpvp(6) - 8*tmpvp(7) + 8*tmpvp(8) - tmpvp(9))/(12*dz);

    % 压强对R和Z的导数
    drdr = (tmpvr(1) - 8*tmpvr(2) + 8*tmpvr(4) - tmpvr(5))/(12*dr);
    drdz = (tmpvr(6) - 8*tmpvr(7) + 8*tmpvr(8) - tmpvr(9))/(12*dz);

% 通量管体积（参考点附近的9点的通量体积）
    tmpvv = 0*tmpvp;
    for k2 = 1:length(tmpvp)
        tmpp  = tmpvp(k2);
        mpfl  = contourc(R(:,1)',Z(1,:),psiout',[tmpp tmpp]);
        % - 保留D线圈的通量(一般来讲就两个)
        indtr = find( mpfl(1,:) < 0 );
        % 方法1 判断那个区域更长 ； 返回位置
        if ( indtr == 1 )
            tmpfl = mpfl(:,2:end);
        else %包括多个分界面取长的
            tmpind = diff([indtr,length(mpfl)]);
            [~,maxtmpind] = max(tmpind);
            tmpfl = mpfl(:,maxtmpind+1:mpfl(2,maxtmpind));
        end

        mflbr  = interp2(R',Z',BBR',tmpfl(1,:),tmpfl(2,:));%插值计算磁场
        mflbz  = interp2(R',Z',BBZ',tmpfl(1,:),tmpfl(2,:));
        % 闭合环积分 -> 周期边界条件: 1.使用第二类曲线积分；2.积分去除顶部L线圈通量
        % 1/B => B/|B|^2 使用磁场分量积分，不是倒数的积分
        % 前后连接的积分
%         amflbr = 0.5*( mflbr + [mflbr(2:end),mflbr(1)] );% R磁场平均值
%         amflbz = 0.5*( mflbz + [mflbz(2:end),mflbz(1)] );%
%         amflbb = amflbr.^2 + amflbz.^2;
%         mfdr = tmpfl(1,:) - [tmpfl(1,2:end),tmpfl(1,1)]; % R位置
%         mfdz = tmpfl(2,:) - [tmpfl(2,2:end),tmpfl(2,1)];
        
        % 收尾不连接的积分
        amflbr = 0.5*( mflbr(1:end-1) + mflbr(2:end) );
        amflbz = 0.5*( mflbz(1:end-1) + mflbz(2:end) );
        amflbb = amflbr.^2 + amflbz.^2;
        mfdr   = tmpfl(1,2:end) - tmpfl(1,1:end-1);
        mfdz   = tmpfl(2,2:end) - tmpfl(2,1:end-1);
        % - 第二类曲线积分
        tmpvv(k2) = sum( mfdr.*amflbr./amflbb + mfdz.*amflbz./amflbb );
    end
% 通量管体积对R和Z的导数（一阶的四点中心差分）
    dvdr = (tmpvv(1) - 8*tmpvv(2) + 8*tmpvv(4) - tmpvv(5))/(12*dr);
    dvdz = (tmpvv(6) - 8*tmpvv(7) + 8*tmpvv(8) - tmpvv(9))/(12*dz);
      
    % 参考点的V和dVdp
    Dvdp(jr) = dvdr./dpdr + dvdz./dpdz; %dVdp
    Dpdp(jr) = drdr./dpdr + drdz./dpdz; %dperdp
    Vf(jr)   = tmpvv(3);
end

