function [Vftv] = FTV(psiout,BBR,BBZ,R,Z,mpp)
%FTV Flux Tube Volume
%   输入 磁通量psi , R 和 Z BR 和 BZ；返回对应的通量管体积;以及计算体积的目标通量和对应的半径

Vftv  = 0*mpp;                   % 通量管体积

for k2 = 1:length(Vftv)
    tmpp  = mpp(k2);
    mpfl  = contourc(R(:,1)',Z(1,:),psiout',[tmpp tmpp]);
    % - 保留D线圈的通量(一般来讲就两个)
    indtr = find( mpfl(1,:) < 0 );
    % 方法1 判断那个区域更长 ； 返回位置
    if ( length(indtr) == 1 )
        tmpfl = mpfl(:,2:end);
    else %包括多个分界面最长的
        tmpind = diff([indtr,length(mpfl)]);
        maxtmpind = find( tmpind == max(tmpind) );
        tmpfl = mpfl(:,indtr(maxtmpind)+1:mpfl(2,indtr(maxtmpind))+indtr(maxtmpind));
    end

    mflbr  = interp2(R',Z',BBR',tmpfl(1,:),tmpfl(2,:),'linear');%插值计算磁场
    mflbz  = interp2(R',Z',BBZ',tmpfl(1,:),tmpfl(2,:),'linear');
%     max(isnan(mflbr))
%     max(isnan(mflbz))
%     % 闭合环积分 -> 周期边界条件: 1.使用第二类曲线积分；2.积分去除顶部L线圈通量
%     amflbr = 0.5*( mflbr + [mflbr(2:end),mflbr(1)] );% 磁场平均值
%     amflbz = 0.5*( mflbz + [mflbz(2:end),mflbz(1)] );
%     amflbb = amflbr.^2 + amflbz.^2;
%     mfdr = tmpfl(1,:) - [tmpfl(1,2:end),tmpfl(1,1)];
%     mfdz = tmpfl(2,:) - [tmpfl(2,2:end),tmpfl(2,1)];
    % 非闭合积分
    amflbr = 0.5*( mflbr(2:end) + mflbr(1:end-1) );% 磁场平均值
    amflbz = 0.5*( mflbz(2:end) + mflbz(1:end-1) );
    amflbb = amflbr.^2 + amflbz.^2;
    mfdr   = diff(tmpfl(1,:));
    mfdz   = diff(tmpfl(2,:));
    % - 第二列曲线积分
    tvfv = mfdr.*amflbr./amflbb + mfdz.*amflbz./amflbb;
    Vftv(k2) = sum( tvfv );
end
end

