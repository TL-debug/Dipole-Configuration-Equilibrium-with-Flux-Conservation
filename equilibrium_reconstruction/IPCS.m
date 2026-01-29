function [tmpP1,tmpJ1,RHSJ,tmpI1] = IPCS(psiout,p0,Ho,Hi,rh,g,mp,R,Z)
% IPSC Initialization of Pressure, Current and Source term
% 初始化等离子体压强，电流和源项
%  输入 ；输出 
mu0 = 4*pi*10^-7;
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);ds =dr*dz;
r1 = R(1,1);z1 = Z(1,1);NR = size(R,1);NZ = size(R,2);
if ( mp == 1 )
    tmpPSIC = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC(tmpPSIC~=0));mnpsic=min(min(tmpPSIC));
    zh = 0;dnh = floor(0.05/dr)+1;%相邻三个基本单位的加热宽度
    % 加热位置
    rhnu=floor((rh-r1)/dr+2);
    zhnu=floor((zh-z1)/dz+1);
    psih=tmpPSIC(rhnu,zhnu);      % 加热位置；加热位置固定
    dpsim=tmpPSIC(rhnu+dnh,zhnu); % 加热位置右侧通量
    dpsip=tmpPSIC(rhnu-dnh,zhnu); % 加热位置左侧通量
    % ---- 3 计算 过度函数系数 A B C （ 加权 Or 平均 ） kap = 0.4
    %g = 2.56;
    alp = 4*g*(abs(mnpsic/psih)-1);
    % g = 2.54;p0 = 1.9801 最大beta在10附近;p0=1.7020,betam=5；g = 2.9 p0=2909 beta = 1000;g=2.8 p0=2.5939e4 beta=100
    % g = 1.0;  p0 = 2455 ; beta = 14w;
    % g = 1.2;  p0 = 2200 ; beta = 170.576;
    % g = 1.4;  p0 = 2125 ; beta = 16.6435;
    % g = 1.6;  p0 = 2250 ; beta = 5.9728
    % g = 1.8;  p0 = 2525 ; beta = 3.6566
    % g = 2.0;  p0 = 2925 ; beta = 3.2561
    % g = 2.2;  p0 = 3425 ; beta = 4.4882
    % g = 2.4;  p0 = 4000;
    % g = 2.6;  p0 = 4700;
    % g = 2.8;  p0 = 5500;
    % g = 3.0;  p0 = 6320;
    % ------------------------------
    MCA=[2*dpsip 1 0;2*dpsim 1 0;dpsip^2 dpsip 1;dpsim^2 dpsim 1];
    MCB=[alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
        4*g*p0/psih*(dpsim/psih)^(4*g-1);
        p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
        p0*(dpsim/psih)^(4*g)];
    %TFC=lsqr(MCA,MCB);% 慢
    %TFC=pinv(MCA)*MCB;% 中
    TFC=MCA\MCB;       % 快
    % ---- 对等离子体区域进行赋值
    indk = find( tmpPSIC <= dpsip );
    tmpP1(indk) = p0*((tmpPSIC(indk)-mnpsic)/(psih-mnpsic)).^alp;
    tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((tmpPSIC(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);
    indk = find( tmpPSIC > dpsip & tmpPSIC <= dpsim );
    tmpP1(indk) = TFC(1)*tmpPSIC(indk).^2+TFC(2)*tmpPSIC(indk)+TFC(3);
    tmpJ1(indk) = -R(indk).*(2*TFC(1)*tmpPSIC(indk)+TFC(2));
    indk = find( tmpPSIC > dpsim );
    tmpP1(indk) = p0*(tmpPSIC(indk)/psih).^(4*g);
    tmpJ1(indk) = -4*g*p0*R(indk)/psih.*(tmpPSIC(indk)/psih).^(4*g-1);
    tmpP1 = reshape(tmpP1,NR,NZ);
    tmpJ1 = reshape(tmpJ1,NR,NZ);
elseif ( mp == 2 )
    % - Third Tpye Pressure  &  Current Profile ( D. T. Garnier 1999 )
    tmpPSIC0 = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC0(tmpPSIC0~=0));mnpsic=min(min(tmpPSIC0));

    tmpP1=0.5*p0*(1-cos(2*pi*(psiout-mxpsic)/(mnpsic-mxpsic)));
    tmpJ1=R.*p0*pi/(mxpsic-mnpsic).*sin(2*pi*(psiout-mxpsic)./(mnpsic-mxpsic));
end
% - 删除非物理区域
tmpP1 = Hi.*Ho.*tmpP1;
tmpJ1 = Hi.*Ho.*tmpJ1;
tmpI1 = sum(sum(tmpJ1))*ds;
RHSJ = mu0*R.*tmpJ1;

end

