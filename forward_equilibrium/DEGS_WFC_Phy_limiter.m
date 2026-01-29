% DEGS _ with flux consersation _ Physics limiter
% This configuration confines the expansion of the LCFS...
% (Last Closed Flux Surface) through physics limiter
% physics limiter is located on the inner side of the vacuum vessel
% ( !Only of forward equilibrium ! ) 
%% ---- magnetic limiter configuration ----- %%
clear;clc;close all;
%  ---- set the computational domain and mesh density settings ---- %
nI = 8;NG = 2^nI + 1;nlay = nI - 1;NR = NG;NZ = NG;% number of grids and multigrid levels
r1 = 0.0;r2 = 4.0;z1 = -1.5;z2 = 2.5;% computational domain
c0 = 2e-7;mu0 = 4*pi*1e-7;           % Vacuum permeability c0 = mu0/2pi
rr = linspace(r1,r2,NR);zz = linspace(z1,z2,NZ);
dr = rr(2) - rr(1);dz = zz(2) - zz(1);dr2 = dr*dr;dz2 = dz*dz;ds = dr*dz;
[R,Z] = ndgrid(rr,zz);

% - coordinate grid at each level of the multigrid method,...
%   along with the corresponding interpolation and restriction matrices - 
[MRI] = MatRI(R,Z,nlay);

% - Green's function between grids with boundary grids
[GrE,indnb,nb] = Gree(rr,zz);%Columns are grid points; rows are boundary points
nbr = indnb(1,:);nbz = indnb(2,:);

% ---- L and D coils ----
% - Coil position and number of turns -
CR=[0.85 0.53];CZ=[2.1 0];LNR = 60;LNZ = 30;DNR = 40;DNZ = 40;
% - Dimensions of the L-coil - 
LW=0.3;LH=0.15;
% - Dimensions of the D-coil - 
DW = 0.16;DH = 0.16;
%  Current of Coils 
CJ = [4.7e5 4.77e6];
% Dipole coil and Levitated coil
VCRD = linspace(CR(2) - DW/2,CR(2) + DW/2,DNR);
VCZD = linspace(CZ(2) - DH/2,CZ(2) + DH/2,DNZ);
VCRL = linspace(CR(1) - LW/2,CR(1) + LW/2,LNR);
VCZL = linspace(CZ(1) - LH/2,CZ(1) + LH/2,LNZ);
[CRD,CZD] = ndgrid(VCRD,VCZD);
[CRL,CZL] = ndgrid(VCRL,VCZL);

% - assemble coils  （ L_Coil D_Coil ）
Coils = [[{CRL},{CZL},{CJ(1)}];[{CRD},{CZD},{CJ(2)}]];
Coils0 = Coils;

% ---- Mutual Inductance Matrix and Self-Inductance ----
[mDL] = FMI(CRD,CZD,CRL,CZL);  % The mutual inductance between D-Coil  and L-Coil
[InD] = FSI(CRD,CZD);          % D-coil self-inductance
[ncd,MuT] = MGC(R,Z,CRD,CZD);  % The mutual inductance between D-Coil and grids:Rows represent coils; columns represent the mesh

% ---- Interpolate the coil currents onto the mesh ----
[RHSC] = p2gv(Coils,rr,zz);
% - compute boundary conditions
psibc = -mu0*GrE*reshape(RHSC,[],1);
% - source term of the G-S
RHSC   = mu0*R.*RHSC/ds;

% ---- initial guess magnetic flux ----
[PSIC] = MCVT(R,Z,Coils);

%% ---- Plasma region and pressure profile initialization
% - Initial magnetic flux through the D-coil
tID  = CJ(2);tIL  = CJ(1);
IT = tID*tIL;
cp0 = InD*CJ(2)/(DNR*DNZ) + mDL*CJ(1)/(LNR*LNZ);
rDC = 0.15;RDC = 0.55; % major and minor radius of D-oil
RSOL = 3.8;            % location limiter (m)

omeo = 0.80;           % SOR acceleration factor
nt   = 200;            % Maximum number of iterations for inner loop
tolo = 10^-9;          % Threshold condition for outer loop
psiout = PSIC;         % guess flux
rh = 1.5;              % location of peak pressure
pp = 0.9e4;dp = 1000;p0 = pp; % peak pressure
tmpi0 = 0.1e6;                % Target plasma current
dcur = 1;                     % Relative error of current
% - 
dpso = 10*tolo;       % Error of the boundary                                            
dcp  = 1;             % Error of the flux through the D-coil
ito  = 1;             % outer loop calculation

% ---- Memory
% - 3D matrixs (pressure, current density, flux)
MPRE = [];
MCUR = [];
MPSI = [];
% - vectors ( Coil currents & peak pressure & plasma beta % flux through D-coil)
IR   = [];
VP0  = [];
VBet = [];
LAM  = [];
% - 2D matrixs ( flux and current density )
% Convert the flux into a column and save it as a 2D matrix, checking... 
% whether the current density and flux change abruptly during the iterations
ncu  = 10; % Record every 10 inner-loop iterations.
VCI  = []; % current density
VPS  = []; % flux
% - 
IDLP = []; % Coil and plasma currents

% ---- dentify interfaces and initialize pressure distribution -----
tmpPSIC = psiout;
[Mo1,Mi1,Ho,Hi,~] = SIPL(tmpPSIC,RSOL,RDC,rDC,R,Z);
Mi=Mi1;Mo=Mo1;% Initial inner and outer plasma edges

% - pressure profiles
mp = 3; % Select pressure distribution
if ( mp == 1 )% Davis "Pressure profiles of plasmas confined in the field of a dipole magnet"
    tmpPSIC = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC(tmpPSIC~=0));mnpsic=min(min(tmpPSIC));
    zh = 0;dnh = floor(0.07/dr)+1;% Width of transition area at middle plane
    rhnu=floor((rh-r1)/dr+2);
    zhnu=floor((zh-z1)/dz+1);
    psih=tmpPSIC(rhnu,zhnu);      % flux on peak pressure
    dpsim=tmpPSIC(rhnu+dnh,zhnu); % Flux on right of the transition area
    dpsip=tmpPSIC(rhnu-dnh,zhnu); % Flux on left of the transition area
    % - Calculate coefficients A B C -
    g = 1.9;alp = 4*g*(abs(mnpsic/psih)-1);
    MCA=[2*dpsip 1 0;2*dpsim 1 0;dpsip^2 dpsip 1;dpsim^2 dpsim 1];
    MCB=[alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
         4*g*p0/psih*(dpsim/psih)^(4*g-1);
         p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
         p0*(dpsim/psih)^(4*g)];
    TFC=MCA\MCB;
    % - pressure ad current density on plasma region -
    indk = find( tmpPSIC <= dpsip );
    tmpP1(indk) = p0*((tmpPSIC(indk)-mnpsic)/(psih-mnpsic)).^(alp);
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

elseif( mp == 3 )
    tmpPSIC = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC(tmpPSIC~=0));mnpsic=min(min(tmpPSIC));
    zh = 0;dnh = floor(0.1/dr)+1;
    rhnu=floor((rh-r1)/dr+2);
    zhnu=floor((zh-z1)/dz+1);
    psih=tmpPSIC(rhnu,zhnu);
    dpsim=tmpPSIC(rhnu+dnh,zhnu);
    dpsip=tmpPSIC(rhnu-dnh,zhnu);
    % - Calculate coefficients A B C D- 
    g = 2.0;alp = 4*g*(abs(mnpsic/psih)-1);
    MCL = [ dpsip^3  dpsip^2 dpsip 1;
            dpsim^3  dpsim^2 dpsim 1;
           3*dpsip^2 2*dpsip   1   0;
           3*dpsim^2 2*dpsim   1   0];
    MCR = [ p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
        p0*(dpsim/psih)^(4*g);
        alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
        4*g*p0/psih*(dpsim/psih)^(4*g-1)];
    TFC = MCL\MCR;
     % - pressure ad current density on plasma region -
    indk = find( tmpPSIC <= dpsip ); % peak pressure inner region
    tmpP1(indk) = p0*((tmpPSIC(indk)-mnpsic)/(psih-mnpsic)).^alp;
    tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((tmpPSIC(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);
    indk = find( tmpPSIC > dpsip & tmpPSIC <= dpsim );% transition zone
    tmpP1(indk) = TFC(1)*tmpPSIC(indk).^3 + TFC(2)*tmpPSIC(indk).^2 + TFC(3)*tmpPSIC(indk) + TFC(4);
    tmpJ1(indk) = -R(indk).*(3*TFC(1)*tmpPSIC(indk).^2 + 2*TFC(2)*tmpPSIC(indk) + TFC(3));
    indk = find( tmpPSIC > dpsim );
    tmpP1(indk) = p0*(tmpPSIC(indk)/psih).^(4*g);
    tmpJ1(indk) = -4*g*p0*R(indk)/psih.*(tmpPSIC(indk)/psih).^(4*g-1);
    tmpP1 = reshape(tmpP1,NR,NZ);
    tmpJ1 = reshape(tmpJ1,NR,NZ);
end
% - identified physical regions
tmpP1 = Hi.*Ho.*tmpP1;
tmpJ1 = Hi.*Ho.*tmpJ1;
tmpI1 = sum(sum(tmpJ1))*ds
RHSJ = mu0*R.*tmpJ1;

% - compute boundary condition (plasma) 
psibp = -mu0*ds*GrE*reshape(tmpJ1,[],1);

% - SOR algorithm register
RHSC1  = RHSC;RHSC2 = RHSC;   % soure term of coils
psibc1 = psibc;psibc2 = psibc;% boundary condition of coils
psibp1 = psibp;psibp2 = psibp;% boundary condition of plasma current
%% ---- major loop ----- 
while( dpso >= tolo || dcp >= 1e-3 || dcur >= 1e-3 || dpsi >= tolo)% boundary flux; flux through D-coil; target current; internal flux 
    
    % - update soure terms
    RHS = RHSC + RHSJ;
    
    % update boundary condition ( plasma + coils )
    psib = psibc + psibp;
    for j = 1:nb
        psiout(nbr(j),nbz(j)) = psib(j);
    end
    
    % - inner loop initialize
    iti  = 1;         % count
    toli = 10^-10;    % threshold
    dpsi = 10*toli;   % error
    omei = omeo;      % SOR acc factor
    NIT  = 4;         % iteration number
    %
    psi0 = psiout;psit = psiout;
    while ( iti <= nt && dpsi >= toli )
        % ------- SOR --------
        [psit,RES] = GSSOR2(psi0,RHS,R,Z,psib,nbr,nbz,NIT,omei);
        dpsi = max(max(abs(psit(2:end-1,2:end-1) - psi0(2:end-1,2:end-1))));
        % ------ MG Acc -------
        [psi0] = MGSOR2(psit,RES,MRI,NIT,omei);
        % ---- Separatrix -----
        [Mo,Ho,~] = SIPLO(psi0,RSOL,R,Z,Hi);
        % ---- Pressure profile ----
        if ( mp == 1 )
            tmpPSIC0 = Hi.*Ho.*psi0;
            mxpsic = max(tmpPSIC0(tmpPSIC0~=0));mnpsic = min(min(tmpPSIC0));
            psih=psi0(rhnu,zhnu);      
            dpsim=psi0(rhnu+dnh,zhnu); 
            dpsip=psi0(rhnu-dnh,zhnu); 
            alp = 4*g*(abs(mnpsic/psih)-1);
            MCA=[2*dpsip 1 0;2*dpsim 1 0;dpsip^2 dpsip 1;dpsim^2 dpsim 1];
            MCB=[alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
                4*g*p0/psih*(dpsim/psih)^(4*g-1);
                p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
                p0*(dpsim/psih)^(4*g)];
            TFC=MCA\MCB;
            indk = find( psi0 >= mnpsic & psi0 <= dpsip );
            tmpP1(indk) = p0*((psi0(indk)-mnpsic)/(psih-mnpsic)).^alp;
            tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((psi0(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);
            indk = find( psi0 > dpsip & psi0 <= dpsim );
            tmpP1(indk) = TFC(1)*psi0(indk).^2+TFC(2)*psi0(indk)+TFC(3);
            tmpJ1(indk) = -R(indk).*(2*TFC(1)*psi0(indk)+TFC(2));
            indk = find( psi0 > dpsim );
            tmpP1(indk) = p0*(psi0(indk)/psih).^(4*g);
            tmpJ1(indk) = -4*g*p0*R(indk)/psih.*(psi0(indk)/psih).^(4*g-1);
        elseif ( mp == 2 )
            tmpPSIC0 = Hi.*Ho.*psi0;
            mxpsic = max(tmpPSIC0(tmpPSIC0~=0));
            mnpsic = min(min(tmpPSIC0));
            tmpP1 = 0.5*p0*(1-cos(2*pi*(psi0-mxpsic)/(mnpsic-mxpsic)));
            tmpJ1 = R.*p0*pi/(mxpsic-mnpsic).*sin(2*pi*(psi0-mxpsic)./(mnpsic-mxpsic));
        elseif ( mp == 3 )
            tmpPSIC0 = Hi.*Ho.*psi0;
            mxpsic = max(tmpPSIC0(tmpPSIC0~=0));mnpsic = min(min(tmpPSIC0));
            psih=psi0(rhnu,zhnu);      
            dpsim=psi0(rhnu+dnh,zhnu); 
            dpsip=psi0(rhnu-dnh,zhnu); 
            alp = 4*g*(abs(mnpsic/psih)-1);
            MCL = [ dpsip^3  dpsip^2 dpsip 1;
                dpsim^3  dpsim^2 dpsim 1;
                3*dpsip^2 2*dpsip   1   0;
                3*dpsim^2 2*dpsim   1   0];
            MCR = [ p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
                p0*(dpsim/psih)^(4*g);
                alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
                4*g*p0/psih*(dpsim/psih)^(4*g-1)];
            TFC = MCL\MCR;
            indk = find( psi0 <= dpsip ); 
            tmpP1(indk) = p0*((psi0(indk)-mnpsic)/(psih-mnpsic)).^alp;
            tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((psi0(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);
            indk = find( psi0 > dpsip & psi0 <= dpsim );
            tmpP1(indk) = TFC(1)*psi0(indk).^3 + TFC(2)*psi0(indk).^2 + TFC(3)*psi0(indk) + TFC(4);
            tmpJ1(indk) = -R(indk).*(3*TFC(1)*psi0(indk).^2 + 2*TFC(2)*psi0(indk) + TFC(3));
            indk = find( psi0 > dpsim );
            tmpP1(indk) = p0*(psi0(indk)/psih).^(4*g);
            tmpJ1(indk) = -4*g*p0*R(indk)/psih.*(psi0(indk)/psih).^(4*g-1);
        end
        
        % ---- Source term ----
        tmpP1 = tmpP1.*Hi.*Ho;
        tmpJ1 = tmpJ1.*Hi.*Ho;
        tmpI1 = sum(sum(tmpJ1))*ds;
        RHSJ  = mu0*R.*tmpJ1;
        
        RHS = RHSC + RHSJ;
        iti = iti + NIT*nlay; %     

    end

    psiout = psi0;% update psi

    % - Record the current and plasma pressure
    IDLP = [IDLP,[tID;tIL;tmpI1]];
    VP0  = [VP0,p0];
    
    % - Update the coil current based on conservation conditions
    pdm = ds*sum(MuT*reshape(tmpJ1,[],1));
    tcp = tID*InD/(DNR*DNZ)  +  tIL*mDL/(LNR*LNZ)  +  pdm;
    dcp = abs(tcp - cp0)/tcp; % error
    lam = cp0/tcp;            % prop factor 
    tID = lam*tID;            % update current of L-coil
    tIL = IT/tID;             % force balance between D-coil gravity and magnetic attraction
    LAM  = [LAM,lam];

    % - compute boundary condition (plasma current)
    psibp2 = psibp1;
    psibp1 = -mu0*ds*GrE*reshape(tmpJ1,[],1);
    psibp = omeo*psibp1 + (1 - omeo)*psibp2;
    
	% - record and update coil currents
    Coils(1:2,3) = [{tIL};{tID}];
    RHSC2   = RHSC1;
    [RHSC1] = p2gv(Coils,rr,zz);
    RHSC    = omeo*RHSC1 + (1 - omeo)*RHSC2; % 适配源在边界值
    % - compute boundary condition (coils current)
    psibc = -mu0*GrE*reshape(RHSC,[],1);
    RHSC   = mu0*R.*RHSC/ds;

    % - error of boundary magnetic flux
    dpso = max( max(abs(psibp2 - psibp1)) , max(abs(psibc2 - psibc1)) );
    
    if ( ito > 3 )% Constrain the maximum allowable pressure
        p0s = tmpi0/IDLP(3,ito)*VP0(ito);
        if ( p0s > 2.0e4 | p0s <= 0 )
            p0s = 2.0e4;
        end
        p0 = (1-omeo)*VP0(ito) + omeo*p0s;
    end

    dcur = abs(IDLP(3,ito) - tmpi0)/abs(tmpi0);% error between the plasma and target currents

    if ( mod(ito,10)  == 0 )
        disp(['outer loop error = ',num2str(dpso)]);
        disp(['p0 = ',num2str(p0)]);
        disp(['plasma current = ',num2str(tmpI1)]);
    end

    ito = ito + 1;
end
%% ---- Diagnostics ----

subplot(2,2,1)% plasma pressure cloud map

[~,ch2]=contour(R,Z,tmpP1,50,'linewidth',1.5);hold on;axis([r1 r2 z1 z2]);colormap jet;
plot(Mi(1,:),Mi(2,:),'m--','linewidth',2);
spo=plot(Mo(1,2:end-1),Mo(2,2:end-1),'m--','linewidth',2);
plot(CRD,CZD,'-r');plot(CRD',CZD','-r');
plot(CRL,CZL,'-r');plot(CRL',CZL','-r');

set(gca,'TickLabelInterpreter','latex',"FontSize",18);
title('Pressure','Interpreter','latex',"FontSize",18);
%xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",18);
ylabel('$Z\ (\rm m)$','Interpreter','latex',"FontSize",18);
Hh1=legend([ch2 spo],'P','Separatrix');
set(Hh1,'Orientation','vertical','Box','off')
set(Hh1,'Interpreter','latex');
set(Hh1,'FontName','Times New Roman','FontSize',16);


subplot(2,2,2)% plasma current density cloud map
contour(R,Z,tmpJ1,50,'linewidth',1.5);hold on;
axis([r1 r2 z1 z2]);colormap jet;
plot(Mi(1,2:end-1),Mi(2,2:end-1),'m--','linewidth',2);
sho=plot(Mo(1,2:end-1),Mo(2,2:end-1),'m--','linewidth',2);
set(gca,'TickLabelInterpreter','latex',"FontSize",20);
title('Current','Interpreter','latex',"FontSize",20);
%xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",20);
%ylabel('$Z\ (\rm m)$','Interpreter','latex',"FontSize",20);


subplot(2,2,3) % plasma pressure profile at mid-plane
FHL=0;
zlt=floor((FHL-z1)/dz);
indMo1=find(Mo(2,:)<0,1);indMo2=find(Mo(2,:)<0,1,'last');% outer separatrix
tmpRo1=Mo(1,indMo1-2:indMo1+1);tmpZo1=Mo(2,indMo1-2:indMo1+1);
tmpRo2=Mo(1,indMo2-1:indMo2+2);tmpZo2=Mo(2,indMo2-1:indMo2+2);
RMo1=interp1(tmpZo1,tmpRo1,FHL);
RMo2=interp1(tmpZo2,tmpRo2,FHL);
indMi1=find(Mi(2,:)<0,1);indMi2=find(Mi(2,:)<0,1,'last');% inner separatrix
tmpRi1=Mi(1,indMi1-2:indMi1+1);tmpZi1=Mi(2,indMi1-2:indMi1+1);
tmpRi2=Mi(1,indMi2-1:indMi2+2);tmpZi2=Mi(2,indMi2-1:indMi2+2);
RMi1=interp1(tmpZi1,tmpRi1,FHL);
RMi2=interp1(tmpZi2,tmpRi2,FHL);
pp1=tmpP1(1:end,zlt);
zsp1=linspace(min(pp1)/1000,max(pp1)/1000,50);
rsp1=RMi1*ones(size(zsp1));
rsp2=RMi2*ones(size(zsp1));
rsp3=RMo1*ones(size(zsp1));
rsp4=RMo2*ones(size(zsp1));
sho=plot(rsp1,zsp1,'m--','linewidth',2);hold on;
plot(rsp2,zsp1,'m--','linewidth',2);
plot(rsp3,zsp1,'m--','linewidth',2);
plot(rsp4,zsp1,'m--','linewidth',2);
ch1=plot(rr(1:end),pp1/1000,'b','linewidth',1.5);hold off;
ax = gca;
xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",20);
ylabel('$P\ (\rm kPa)$','Interpreter','latex',"FontSize",20);
ax.TickLabelInterpreter = 'latex';
set(gca,'FontName','Times','FontSize',18);
set(gca,'TickLabelInterpreter','latex');

Hh1=legend([ch1 sho],'P','Sepa');
set(Hh1,'Orientation','vertical','Box','off')
set(Hh1,'Interpreter','latex');
set(Hh1,'FontName','Times New Roman','FontSize',16);


subplot(2,2,4)% plasma current density  profile at mid-plane
jj1=tmpJ1(1:end,zlt);
zsp1=linspace(min(jj1)/1000,max(jj1)/1000,50);
rsp1=RMi1*ones(size(zsp1));
rsp2=RMi2*ones(size(zsp1));
rsp3=RMo1*ones(size(zsp1));
rsp4=RMo2*ones(size(zsp1));
sho=plot(rsp1,zsp1,'m--','linewidth',2);hold on;axis tight
plot(rsp2,zsp1,'m--','linewidth',2);
plot(rsp3,zsp1,'m--','linewidth',2);
plot(rsp4,zsp1,'m--','linewidth',2);
ch1=plot(rr(1:end),jj1/1000,'b','linewidth',1.5);hold off;
ax = gca;
xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",20);
ylabel('$J\ (\rm kA/m^2)$','Interpreter','latex',"FontSize",20);
ax.TickLabelInterpreter = 'latex';
set(gca,'FontName','Times','FontSize',18);
set(gca,'TickLabelInterpreter','latex');

Hh1=legend([ch1 sho],'J','Sepa');
set(Hh1,'Orientation','vertical','Box','off')
set(Hh1,'Interpreter','latex');
set(Hh1,'FontName','Times New Roman','FontSize',16);
%% ---- beta cloud map ----

ttmpPsi = psiout;
tmpBR = zeros(size(R));tmpBZ = tmpBR;tmpBB = tmpBR;
for i = 1:NR
    for j = 1:NZ
        if ( j>=3 && j<=NZ-2 && i>=3 && i<=NR-2 )%4点差分
            tmpBR(i,j) =-(1/12*ttmpPsi(i,j+2)-2/3*ttmpPsi(i,j+1)+2/3*ttmpPsi(i,j-1)-...
                          1/12*ttmpPsi(i,j-2))/(dz*R(i,j));
            tmpBZ(i,j) = (1/12*ttmpPsi(i+2,j)-2/3*ttmpPsi(i+1,j)+2/3*ttmpPsi(i-1,j)-...
                         1/12*ttmpPsi(i-2,j))/(dr*R(i,j));
        elseif ( j>=2 && j<=NZ-1 && i>=2 && i<=NR-1 )
            tmpBR(i,j)=-(ttmpPsi(i,j+1)-ttmpPsi(i,j-1))/(2*dz*R(i,j));
            tmpBZ(i,j)= (ttmpPsi(i+1,j)-ttmpPsi(i-1,j))/(2*dr*R(i,j)); 
        end
    end
end

BBG = sqrt(tmpBR.^2 + tmpBZ.^2);
BBG1=BBG.*Hi.*Ho;
BZ=tmpBZ.*Hi.*Ho;
BR=tmpBR.*Hi.*Ho;

% ---- pressure and magnetic filed at mid-plane
% zhnu=floor((0-z1)/dz+1);%mid-plane
%plot(rr(100:end),tmpP1(100:end,zhnu),'-b');hold on;
%plot(rr(100:end),tmpP1(100:end,zhnu)+BZ(100:end,zhnu).^2/(2*mu0),'--m')
% plot(rr(100:end),abs(BZ(100:end,zhnu)),'-r')

beta=2*mu0*tmpP1./BBG1.^2;
beta(isnan(beta))=0;
figure(3)
mesh(R,Z,beta)
abeta=sum(sum(beta))/sum(sum(beta~=0));

maxbeta = max(max(beta))

mesh(R,Z,beta);axis([r1 r2 z1 z2 min(min(beta)) max(max(beta))])

%% ---- coil and plasma currents vary with iteration
h1=plot((IDLP(1,:)-CJ(2))/CJ(2),'-r','linewidth',2);hold on;
h2=plot((IDLP(2,:)-CJ(1))/CJ(2),'--b','linewidth',2);
h3=plot((IDLP(3,:)-IDLP(3,1))/CJ(2),'-.g','linewidth',2);
plot(1:length(IDLP(1,:)),zeros(size(IDLP(1,:))),':m','linewidth',2);
axis tight;grid on; box on ;grid minor;
L1=legend([h1 h2 h3],{'$I_D$', '$I_L$','$I_p$'},"Interpreter","latex","FontSize",18,'FontName','times');
set(gca,'TickLabelInterpreter','latex','fontsize',18);
set (L1,'box','off')
xlabel('Iteration Number','fontname','times')
ylabel('Relative Current Change Rate','fontname','times')

%% ---- Output the data required for equilibrium reconstruction(Not applicable to equilibrium inverse solutions)
% % - Select fixed points along the plasma boundary
% plot(Mo(1,2:end),Mo(2,2:end),'.m','linewidth',2);axis([r1 r2 z1 z2]);hold on;
% clear GrE;clear MuT;
% [~,beta1,~] = GradB(psiout,tmpP1,R,Z);
% % - identify separatrix
% indsp = find(Mo(1,:)<0);
% indsp = [indsp,length(Mo)];[~,indspo] = max(diff(indsp));
% indmo = indsp(indspo) + 1:indsp( indspo + 1 ) - 1;
% % - compute magnetic field along separatrix
% tmpsep = [Mo(1,indmo);Mo(2,indmo)];
% [FBR,FBZ,DFBR,DFBZ] = GradFB(tmpsep,psiout,R,Z,dr,dz);
% [~,indmbb] = min(FBR.^2 + FBZ.^2);
% % - determine fixed points (select fixed points at equal intervals along separatrix)
% np = 20;% number of fixed point
% indtmp  = indmo(1:floor(length(indmo)/np):end);
% indtmp1 = indtmp+indmbb-indtmp(1);             
% tmptmpind = find(indtmp1 > indtmp(end));       
% indtmp1(tmptmpind) = indtmp1(tmptmpind) - indtmp(end);
% % - first and second rows correspond to the R and Z of fixed points
% Sep = [tmpsep(1,indtmp1(1:end-1));tmpsep(2,indtmp1(1:end-1))];
% 
% % Calculate flux and transverse and longitudinal magnetic field at fixed points
% psisep = interp2(R',Z',psiout',Sep(1,:),Sep(2,:));
% [FBR,FBZ,DFBR,DFBZ] = GradFB(Sep,psiout,R,Z,dr,dz); 
% Sep = [Sep;psisep;FBR;FBZ;DFBR;DFBZ]; % add fixed points information
% 
% % plot fixed points
% plot(Sep(1,1),Sep(2,1),'xg','markersize',15,'linewidth',3);%X点位置
% plot(Sep(1,2:end),Sep(2,2:end),'xb','markersize',15,'linewidth',3);%其他固定点
% mesh(R,Z,psiout);hold on;
% plot3(Sep(1,:),Sep(2,:),Sep(3,:),'or','markersize',15,'linewidth',3);
% grid on;
% set(gca,'FontName','Times','FontSize',18);
% set(gca,'TickLabelInterpreter','latex');
% xlabel('$R\ (\rm m)$','Interpreter','latex',"FontSize",18);
% ylabel('$Z\ (\rm m)$','Interpreter','latex',"FontSize",18);
% zlabel('$\psi_{plasma}\ (\rm m)$','Interpreter','latex',"FontSize",18);
% title({'\psi & and Fixed plasma boundary position '})
% 
% % ---- output data for equilibrium reconstruction
% CPLA =  cat(3,psiout,tmpP1,tmpJ1,beta1);%flux,pressure,current density,beta
% save(['sep_points_phy_limiter',num2str(np),'.mat'],'Sep','Mo','IDLP','Coils0','CPLA');

