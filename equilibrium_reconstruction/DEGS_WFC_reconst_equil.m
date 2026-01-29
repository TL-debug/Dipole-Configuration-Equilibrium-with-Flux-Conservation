% Update:2024/12/1 _ by:liu T
% DEGS_with flux consersation_equilibrium reconstruction
% The solver utilizes the positions and measurement data of magnetic flux...
% measurement coils and can reconstruct the equilibrium by applying constraints...
% on physical quantities (such as plasma current, plasma pressure, beta value, etc.)....
% Additionally, it is also capable of performing equilibrium reconstruction without...
% introducing such physical constraints.
% -------------------------------------------------------------------------
%  Reconstruction results of levitated dipole configuration 
%  Result: ID = 4734834.8811(reconstruction);ID0 = 4734853.0543(free boundary)
%          IL = 473423.1051 (reconstruction);IL0 = 473488.823  (free boundary)
%          Ip = 155617.9503 (reconstruction);Ip0 = 155617.9513 (free boundary)
%          P0 =  18054.7737 (reconstruction);P0  =1.8e4        (free boundary)
% -------------------------------------------------------------------------

%% ---- equilibrium reconstruction ----- %%
clear;clc;close all;
% ---------- Initial flux subprogram ----------- % 
% initial Sys 
nI = 8;NG = 2^nI + 1;nlay = nI - 1;NR = NG;NZ = NG;% number of grids and multigrid levels
r1 = 0.0;r2 = 4.5;z1 = -1.7;z2 = 2.5;% computational domain
c0 = 2e-7;mu0 = 4*pi*1e-7;           % Vacuum permeability c0 = mu0/2pi
rr = linspace(r1,r2,NR);zz = linspace(z1,z2,NZ);
dr = rr(2) - rr(1);dz = zz(2) - zz(1);dr2 = dr*dr;dz2 = dz*dz;ds = dr*dz;
[R,Z] = ndgrid(rr,zz);

% - coordinate grid at each level of the multigrid method,...
%   along with the corresponding interpolation and restriction matrices - 
[MRI] = MatRI(R,Z,nlay);

% - Green's function between grids with boundary grids
[GrE,indnb,nb] = Gree(rr,zz);% Columns are grid points; rows are boundary points
nbr = indnb(1,:);nbz = indnb(2,:);

% ---- L and D coils ----
% - Coil position and number of turns -
CR=[0.85 0.53];CZ=[2.1 0];LNR = 60;LNZ = 30;DNR = 40;DNZ = 40;
% - Dimensions of the L-coil - 
LW=0.3;LH=0.15;
% - Dimensions of the D-coil - 
DW = 0.16;DH = 0.16;
% - Current of Coils 
CJ = [4.7e5 4.77e6];
% - Dipole coil and Levitated coil
VCRD = linspace(CR(2) - DW/2,CR(2) + DW/2,DNR);
VCZD = linspace(CZ(2) - DH/2,CZ(2) + DH/2,DNZ);
VCRL = linspace(CR(1) - LW/2,CR(1) + LW/2,LNR);
VCZL = linspace(CZ(1) - LH/2,CZ(1) + LH/2,LNZ);
[CRD,CZD] = ndgrid(VCRD,VCZD);
[CRL,CZL] = ndgrid(VCRL,VCZL);

% ---- Reading the input data for equilibrium reconstruction ---- %
% % - input data is selected along separatrix - 
% load('input_data_sep_points_low_divertor_20.mat'); 
% load('input_data_sep_points_mag_limiter_20.mat'); 
% % - input data is selected along position of potential dector - 
% load('input_data_det_points_low_divertor_20.mat'); 
% load('input_data_det_points_mag_limiter_20.mat'); 
% - input data from positions of the flux measurement coils of dipole device - 
load('input_data_Junior_dector_12.mat'); 
% -------------------------------------------------------------------------

% ---- set initial coil currents ---- 
IDLP0 = IDLP;
tmpi0 = 1.0*IDLP0(3,end);% 1.2 times the target current
Coils = Coils0;
% - initial currents of auxiliary coil
for jce = 3:size(Coils,1)
    Coils(jce,3)={-1};
end

% ---- Mutual Inductance Matrix and Self-Inductance ----
[mDL] = FMI(CRD,CZD,CRL,CZL);  % the mutual inductance between D-Coil  and L-Coil
[InD] = FSI(CRD,CZD);          % D-coil self-inductance
[ncd,MuT] = MGC(R,Z,CRD,CZD);  % The mutual inductance between D-Coil and grids:Rows represent coils; columns represent grid

% ---- Interpolate the coil currents onto the mesh ----
[RHSCt] = p2gv(Coils,rr,zz);
% - Boundary conditions
psibct  = -mu0*GrE*reshape(RHSCt,[],1);
RHSCt   =  mu0*R.*RHSCt/ds;

% ---- initial guess magnetic flux ----
[PSIC]  = MCVT(R,Z,Coils);

% ---- GrP Green's function between coil and observation point ---- 
[GrP]   = Grep(Coils,Sep);% mu0*GrP*cell2mat(Coils(:,3))

% ---- Select n observation points along the separatrix ----- %
% - identify separatrix
indsp = find(Mo(1,:)<0);
indsp = [indsp,length(Mo)];[~,indspo] = max(diff(indsp));
indmo = indsp(indspo) + 1:indsp( indspo + 1 ) - 1;
% - compute magnetic field along separatrix
tmpsep = [Mo(1,indmo);Mo(2,indmo)];
[FBR,FBZ,DFBR,DFBZ] = GradFB(tmpsep,CPLA(:,:,1),R,Z,dr,dz);
[~,indmbb] = min(min(FBR.^2 + FBZ.^2));
% - determine fixed points (select fixed points at equal intervals along separatrix)
np = 20;% number of fixed point
indtmp    = indmo(1:floor(length(indmo)/np):end);
indtmp1   = indtmp+indmbb-indtmp(1);            
tmptmpind = find(indtmp1 > indtmp(end));       
indtmp1(tmptmpind) = indtmp1(tmptmpind) - indtmp(end);
Sep0 = [tmpsep(1,indtmp1(1:end-1));tmpsep(2,indtmp1(1:end-1))];%
% ----------------------------------------------------------------------- %
%% ---- 初始化
% ---- colormap for the evolution of coil current ---- 
ColM = slanCL(3,1:size(Coils,1)+1);% number of coils and plasma current

% - Initial magnetic flux through the D-coil
cp0 = InD*CJ(2)/(DNR*DNZ) + mDL*CJ(1)/(LNR*LNZ);% Initial total flux

% - parameters setting
rDC = 0.15;RDC = 0.55;   % major and minor radius of dipole coil
rh = 1.5;                % location of peak pressure

% ---- dentify interfaces and initialize pressure distribution -----
psiout  = PSIC;
tmpPSIC = psiout;
[Mo0,Mi0,Ho,Hi,~] = SIPS(tmpPSIC,RDC,rDC,R,Z); 
Mi = Mi0;Mo = Mo0;

% - guess a peak pressure
p0 = max(max(CPLA(:,:,2)))

% - pressure profiles
mp = 3;% Select pressure distribution
if ( mp == 1 )% Davis "Pressure profiles of plasmas confined in the field of a dipole magnet"
    tmpPSIC = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC(tmpPSIC~=0));mnpsic=min(min(tmpPSIC));
    zh = 0;dnh = floor(0.07/dr)+1;% Width of transition area at middle plane
    rhnu=floor((rh-r1)/dr+2);
    zhnu=floor((zh-z1)/dz+1);
    psih=tmpPSIC(rhnu,zhnu);      % flux on peak pressure
    dpsim=tmpPSIC(rhnu+dnh,zhnu); % flux on right of the transition area
    dpsip=tmpPSIC(rhnu-dnh,zhnu); % flux on left of the transition area
    % - Calculate coefficients A B C -
    g = 2.56 ;alp = 4*g*(abs(mnpsic/psih)-1);
    MCA = [2*dpsip 1 0;2*dpsim 1 0;dpsip^2 dpsip 1;dpsim^2 dpsim 1];
    MCB = [alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
        4*g*p0/psih*(dpsim/psih)^(4*g-1);
        p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
        p0*(dpsim/psih)^(4*g)];
    TFC = MCA\MCB;
    % - pressure and current density on plasma region -
    indk = find( tmpPSIC <= dpsip );
    tmpP1(indk) = p0*((tmpPSIC(indk)-mnpsic)/(psih-mnpsic)).^alp;
    tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((tmpPSIC(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);
    indk = find( tmpPSIC > dpsip & tmpPSIC <= dpsim );% 中间区域赋值
    tmpP1(indk) = TFC(1)*tmpPSIC(indk).^2+TFC(2)*tmpPSIC(indk)+TFC(3);
    tmpJ1(indk) = -R(indk).*(2*TFC(1)*tmpPSIC(indk)+TFC(2));
    indk = find( tmpPSIC > dpsim );  % 外侧区域赋值
    tmpP1(indk) = p0*(tmpPSIC(indk)/psih).^(4*g);
    tmpJ1(indk) = -4*g*p0*R(indk)/psih.*(tmpPSIC(indk)/psih).^(4*g-1);
    tmpP1 = reshape(tmpP1,NR,NZ);
    tmpJ1 = reshape(tmpJ1,NR,NZ);
elseif ( mp == 2 )
    % - Third Tpye Pressure  &  Current Profile ( D. T. Garnier 1999 )
    tmpPSIC0 = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC0(tmpPSIC0~=0));mnpsic=min(min(tmpPSIC0));

    tmpP1=0.5*p0*(1-cos(2*pi*(psiout-mxpsic)/(mnpsic-mxpsic)));
    tmpJ1=R.*p0*pi/(mxpsic-mnpsic).*sin(2*pi*(psiout-mxpsic)./(mnpsic-mxpsic));

elseif ( mp == 3 )% 一元三次的过度函数
    tmpPSIC = Hi.*Ho.*psiout;mxpsic=max(tmpPSIC(tmpPSIC~=0));mnpsic=min(min(tmpPSIC));
    zh = 0;dnh = floor(0.1/dr)+1;
    rhnu=floor((rh-r1)/dr+2);
    zhnu=floor((zh-z1)/dz+1);
    psih=tmpPSIC(rhnu,zhnu);
    dpsim=tmpPSIC(rhnu+dnh,zhnu);
    dpsip=tmpPSIC(rhnu-dnh,zhnu);
    % - Calculate coefficients A B C D- 
    g = 2.0 ;alp = 4*g*(abs(mnpsic/psih)-1);
    MCL = [ dpsip^3  dpsip^2 dpsip 1;
        dpsim^3  dpsim^2 dpsim 1;
        3*dpsip^2 2*dpsip   1   0;
        3*dpsim^2 2*dpsim   1   0];
    MCR = [ p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
        p0*(dpsim/psih)^(4*g);
        alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
        4*g*p0/psih*(dpsim/psih)^(4*g-1)];
    TFC = MCL\MCR;
    % - initialize pressure and current density within plasma region -
    indk = find( tmpPSIC <= dpsip ); % 峰值内侧区域
    tmpP1(indk) = p0*((tmpPSIC(indk)-mnpsic)/(psih-mnpsic)).^alp;
    tmpJ1(indk) = -alp*p0/(psih-mnpsic)*R(indk).*((tmpPSIC(indk)-mnpsic)./(psih-mnpsic)).^(alp-1);

    indk = find( tmpPSIC > dpsip & tmpPSIC <= dpsim );% 过渡区域
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
tmpI1 = sum(sum(tmpJ1))*ds;
RHSJ = mu0*R.*tmpJ1;
RHSC = RHSCt;

% - compute boundary condition (plasma) 
psibp = -mu0*ds*GrE*reshape(tmpJ1,[],1);
psibc = psibct';

% ---- solving the minimum value problem ---- %
% - set regularization coefficient - 
% gamma1 = 0.1e-18;
gamma1 = 0.0005e-15;% 偏滤器位形
% gamma1 = 0.01e-15; % 磁限制器位形
% - compute coefficient matrix - 
MA = 2*mu0^2*(GrP'*GrP) + 2*gamma1*eye(size(GrP,2));
MC = [mDL/(LNR*LNZ) InD/(DNR*DNZ) zeros(1,size(Coils,1)-2)];
MB = MC';
MD = zeros(size(MC,1),size(MC,1));
ML = [MA,MB;MC,MD];
% ----------------------------------------------------------------------- %

% - magnetic flux generated by the coil at the observation points - 
psifc  = mu0*GrP*cell2mat(Coils(:,3));

% ---- Initialize the iterative system parameters.
% - inner loop - 
toli = 1e-9; % threshold
nti  = 100;  % maximum iterations
omei = 0.7;  % SOR acceleration factor
NIT  = 4;    % SOR iteration count
iti  = 1;    % count
% - outer loop -
tolo = 1e-9; % threshold
omeo = 0.7;  % SOR acceleration factor
nito = 500;  % maximum iterations
ito  = 1;    % count
dpso = 0.1;  % error of boundary
dcur = 1;    % error of relative plasma current 

% - SOR algorithm register
tpsibc = psibc;
tpsibp = psibp;
tpsifc = psifc;
VCI1 = cell2mat(Coils(:,3));VCI2 = VCI1;VCI3 = VCI1;
PCD1 = tmpJ1;PCD2 = tmpJ1;PCD3 = tmpJ1;

% ---- Memory
% - 2D matrixs ( peak pressure & beta & plasma current density)
VP0  = [];
VBet = [];
VCur = [];
% - 1D vectors
ncu  = 10; % Record every 10 inner-loop iterations.
VCI  = []; % 
VPS  = []; % psi 
ILDP = []; % currents of palsma and coil

%% ---- major loop -----
% - iteration termination criteria include:
%   1) maximum number of cycles
%   2) flux error in the outer loop
%   3) deviation between the plasma current and the target current
%   4) flux error in the inner loop

while ( ito <= nito && dpso >= tolo || dcur >= 0.1 || dpsi >= toli )

    % - update soure terms
    RHS = RHSC + RHSJ;

    % update boundary condition
    psib = psibc + psibp;
    for j = 1:nb
        psiout(nbr(j),nbz(j)) = psib(j);
    end

    % - inner loop initialize
    iti  = 1;     % conut
    dpsi = 0.1;   % error
    % - 
    psi0 = psiout;psit = psiout;

    % ---- inner loop ---- %
    while ( iti <= nti && dpsi >= toli )
        % ------- SOR --------
        [psit,RES] = GSSOR2(psi0,RHS,R,Z,psib,nbr,nbz,NIT,omei);
        dpsi = max(max(abs(psit(2:end-1,2:end-1) - psi0(2:end-1,2:end-1))));
        % ------ MG Acc -------
        [psi0] = MGSOR2(psit,RES,MRI,NIT,omei);
        % ---- Separatrix -----
        [Mo,Ho,~] = SIPSO(psi0,R,Z,Hi);
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
            % ---- 
            alp = 4*g*(abs(mnpsic/psih)-1);
            MCL = [ dpsip^3  dpsip^2 dpsip 1;
                    dpsim^3  dpsim^2 dpsim 1;
                    3*dpsip^2 2*dpsip   1   0;
                    3*dpsim^2 2*dpsim   1   0];
            MCR = [ p0*((dpsip-mnpsic)/(psih-mnpsic))^alp;
                    p0*(dpsim/psih)^(4*g);
                    alp*p0/(psih-mnpsic)*((dpsip-mnpsic)/(psih-mnpsic))^(alp-1);
                    4*g*p0/psih*(dpsim/psih)^(4*g-1)];
            TFC = MCL\MCR;% 
            % ---- 
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
        iti = iti + NIT*nlay;
    end
    % - update psi
    psiout = psi0; 

    % - Record the coil currents 
    ILDP = [ILDP,[VCI1;tmpI1]];

    % - update boundary condition ( plasma )
    PCD3 = PCD2;PCD2 = tmpJ1;PCD1 = omeo*PCD2 + (1 - omeo)*PCD3;
    RHSJ = mu0*R.*PCD1;
    tpsibp = -mu0*ds*GrE*reshape(PCD1,[],1);

    % - update coil currents
    VCI3 = VCI2;
    % - interpolation flux at observation points（follow updates VHM）
    psift = interp2(R',Z',psiout',Sep(1,:),Sep(2,:),'cubic'); 
    % - flux generated by the coils at observation points
    psifc = mu0*GrP*VCI3;
    psifp = psift' - psifc;

    % - Solve for coil currents via Lagrange multipliers - 
    MR = [2*mu0*GrP'*(Sep(3,:)' - psifp);cp0 - sum(ds*MuT*reshape(PCD1,[],1))];
    tVCI2 = ML\MR;
    VCI2  = tVCI2(1:size(Coils,1));
    VCI1 = omeo*VCI2 + (1 - omeo)*VCI3;

    % - record and update coil currents
    Coils(:,3) = num2cell(VCI1);
    [RHSC] = p2gv(Coils,rr,zz);
    tpsibc = -mu0*GrE*reshape(RHSC,[],1);
    RHSC = mu0*R.*RHSC/ds;

    % - error of boundary magnetic flux 
    dpso = max( [max( abs( tpsibc - psibc ) ),max( abs( tpsibp - psibp ) )] );

    % - save peak pressure and plasma current
    VP0  = [VP0,p0];
    VCur = [VCur,tmpI1];

    % - Constrain the maximum allowable pressure - 
    if ( ito > 3)
        p0 = tmpi0/VCur(end)*VP0(ito);
        if ( p0 > 2e4 || p0 < 0 )
            p0 = 1.99e4;
        end
    end
    % - error between the plasma and target currents
    dcur = abs(VCur(ito) - tmpi0)/tmpi0;

    % ---- visualization of parameters evolution during iterations ----  %
    % - update every 5 iterations
    if ( mod(ito,5) == 0 || ito == 1 )
        subplot(2,2,1)
        % - flux contours
        [~,ch1] = contour(R,Z,psiout,30,'linewidth',2);hold on;
        % - position of flux measurement coils 
        h1 = plot(Sep(1,:),Sep(2,:),'.k','markersize',30); 
        % - separatrix
        indmo = find(Mo(1,:)<0);indmo = [indmo,length(Mo(1,:))];
        for kp = 1:length(indmo)-1
            h2 = plot(Mo(1,indmo(kp)+1:indmo(kp+1)-1),Mo(2,indmo(kp)+1:indmo(kp+1)-1),'--m','linewidth',3);
        end
        % - initial separatrix
        indmo = find(Mo0(1,:)<0);indmo = [indmo,length(Mo0(1,:))];
        for kp = 1:length(indmo)-1
            h3 = plot(Mo0(1,indmo(kp)+1:indmo(kp+1)-1),Mo0(2,indmo(kp)+1:indmo(kp+1)-1),':g','linewidth',3);
        end
        plot(Mi(1,2:end),Mi(2,2:end),'--m','linewidth',3);
        % - target 
        h4 = plot(Sep0(1,:),Sep0(2,:),'.c','markersize',20);

        % - plot coils 
        cola = 0;% adjust the gray level of the coil color
        for jcc = 1:size(Coils,1)
            tCRD = cell2mat(Coils(jcc,1));tCZD = cell2mat(Coils(jcc,2));

            plot(tCRD(1,:),tCZD(1,:),'-','linewidth',1.5,'color',cola*[1 1 1]);
            plot(tCRD(end,:),tCZD(end,:),'-','linewidth',1.5,'color',cola*[1 1 1]);
            plot(tCRD(:,1),tCZD(:,1),'-','linewidth',1.5,'color',cola*[1 1 1]);
            plot(tCRD(:,end),tCZD(:,end),'-','linewidth',1.5,'color',cola*[1 1 1]);
            tmpDR1 = linspace(min(min(tCRD)),max(max(tCRD)),20);
            tmpDZ1 = linspace(min(min(tCZD)),max(max(tCZD)),20);
            plot(tmpDR1,tmpDZ1,'-','linewidth',1.5,'color',cola*[1 1 1]);
            plot(tmpDR1,flip(tmpDZ1),'-','linewidth',1.5,'color',cola*[1 1 1]);
        end
        % - legend - 
        Hh2=legend([h1 h4 h2 h3],'$ob-pots$','$tg-pots$','$sep$','$in-sep$');
        set(Hh2,'Orientation','vertical','Box','off','NumColumns',2)
        set(Hh2,'Interpreter','latex');
        set(Hh2,'FontName','Times New Roman','FontSize',10);

        axis([r1 r2 z1 z2]);
        set(gca,'TickLabelInterpreter','latex',"FontSize",16);
        ylabel('$\rm\ Z\ (\rm m)$','Interpreter','latex',"FontSize",16);
        xlabel('$\rm\ R\ (\rm m)$','Interpreter','latex',"FontSize",16);
        
        drawnow;
        hold off;


        % - errors 
        subplot(2,2,2)% evolution of error with iterations

        error1 = [toli*0.05,max( abs( tpsibp - psibp ) )]; % error in the plasma current at the boundary flux
        error2 = [toli*0.05,max( abs( tpsibc - psibc ) )]; % error in the coil currents at the boundary flux
        error3 = [toli*0.05,max(abs(psift - Sep(3,:)))];   % maximum error of the reconstruction flux with respect to the target value at the observation points.
        error4 = [toli*0.05,sum((psift - Sep(3,:)).^2) + gamma1*sum(VCI1.^2)]; % epsilon

        semilogy(1*ones(size(error1)),error1,'linewidth',20,'color',[0 0 1]);hold on;
        semilogy(2*ones(size(error2)),error2,'linewidth',20,'color',[0 1 0]);
        semilogy(3*ones(size(error3)),error3,'linewidth',20,'color',[1 0 0]);
        semilogy(4*ones(size(error3)),error4,'linewidth',20,'color',[0 1 1]);
        semilogy(0.5:0.1:4.5,toli*ones(size(0.5:0.1:4.5)),'--m','linewidth',3);

        set(gca,'TickLabelInterpreter','latex',"FontSize",16);
        xticks([1 2 3 4])
        xticklabels({'$\psi_{bp}$','$\psi_{bc}$','$\psi_{fp}$','$\epsilon$'})
        ylabel('$\rm\ error$','Interpreter','latex',"FontSize",16);

        % - display iterations
        text(2.2,2*max([max(error1),max(error2),max(error3),max(error4)]),...
            ['INum = ',num2str(ito)],'FontSize',12,'FontName','times');

        axis([0.5 4.5 10^-10 3*max([max(error1),max(error2),max(error3),max(error4)])])

        drawnow;
        hold off;

        subplot(2,2,[3,4])% convergence curve of coil currents
        niterd = 10;      % increase the display range every 10 iterations.
        iterx = 1:niterd*(floor(ito/niterd)+1);
        
        plot(iterx,ones(size(iterx)),':m','linewidth',2);hold on;
        % - D and L coils 
        h1 = plot(ILDP(1,:)/ILDP(1,1),'-','linewidth',2,'color',ColM(1,:));hold on;
             plot(ILDP(1,:)/ILDP(1,1),'.','markersize',20,'color',ColM(1,:));
        h2 = plot(ILDP(2,:)/ILDP(2,1),'-','linewidth',2,'color',ColM(2,:));
             plot(ILDP(2,:)/ILDP(2,1),'o','markersize',10,'linewidth',2,'color',ColM(2,:));        
        % - aux coils
        h3 = plot(ILDP(3,:)/ILDP(3,1),'-','linewidth',2,'color',ColM(3,:));
             plot(ILDP(3,:)/ILDP(3,1),'s','markersize',10,'linewidth',2,'color',ColM(3,:));
        % ---- Junior 
%         h4 = plot(ILDP(4,:)/ILDP(4,1),'-','linewidth',2,'color',ColM(4,:));hold on;
%              plot(ILDP(4,:)/ILDP(4,1),'+','markersize',10,'linewidth',2,'color',ColM(4,:));
%         h5 = plot(ILDP(5,:)/ILDP(5,1),'-','linewidth',2,'color',ColM(5,:));hold on;
%              plot(ILDP(5,:)/ILDP(5,1),'diamond','markersize',10,'linewidth',2,'color',ColM(5,:));
%         h6 = plot(ILDP(6,:)/ILDP(6,1),'-','linewidth',2,'color',ColM(6,:));hold on;
%              plot(ILDP(6,:)/ILDP(6,1),'^','markersize',10,'linewidth',2,'color',ColM(6,:));
        % ---- magnetic limiter and divertor
%         h7 = plot(ILDP(7,:)/ILDP(7,1),'-','linewidth',2,'color',ColM(7,:));hold on;
%              plot(ILDP(7,:)/ILDP(7,1),'pentagram','markersize',20,'linewidth',2,'color',ColM(7,:));
%         h8 = plot(ILDP(8,:)/ILDP(8,1),'-','linewidth',2,'color',ColM(8,:));hold on;
%              plot(ILDP(8,:)/ILDP(8,1),'x','markersize',10,'linewidth',2,'color',ColM(8,:));
%         h9 = plot(ILDP(9,:)/ILDP(9,1),'-','linewidth',2,'color',ColM(9,:));hold on;
%              plot(ILDP(9,:)/ILDP(9,1),'v','markersize',10,'linewidth',2,'color',ColM(9,:));
%         h10 = plot(ILDP(10,:)/ILDP(10,1),'-','linewidth',2,'color',ColM(10,:));
%              plot(ILDP(10,:)/ILDP(10,1),'s','markersize',10,'linewidth',2,'color',ColM(10,:));
        % ---- Up-down symmetric divertor 
        % - set axis and legend 
%         Hh1=legend([h1 h2 h3 h4 h5 h6],'$I_{L}$','$I_{D}$','$I_{o1}$','$I_{o2}$','$I_{o3}$','$I_{P}$');
        Hh1=legend([h1 h2 h3],'$I_{L}$','$I_{D}$','$I_{P}$');
        set(Hh1,'Orientation','vertical','Box','off')
        set(Hh1,'Interpreter','latex');
        set(Hh1,'FontName','Times New Roman','FontSize',16);

        set(gca,'TickLabelInterpreter','latex',"FontSize",16);
        ylabel('$\rm\ Relative\ current$','Interpreter','latex',"FontSize",16);
        xlabel('$\rm\ iteration$','Interpreter','latex',"FontSize",16);

        axis([1 max(iterx) 0.95 1.25]);
        drawnow;
        hold off;
    end


    % - update boundary condition (plasma + coils)
    psibp = tpsibp; 
    psibc = tpsibc; 

    ito = ito + 1;
end

% - difference between equilibrium reconstruction and input parameters -
disp(['outer loop error = ',num2str(dpso)]);
disp(['peak pressure = ',num2str(p0)]);
disp(['ID = ',num2str(VCI2(2)),'; ID0 = ',num2str(IDLP0(1,end))]);
disp(['IL = ',num2str(VCI2(1)),'; IL0 = ',num2str(IDLP0(2,end))]);
% disp(['Io1 = ',num2str(VCI2(3)),'; Io10 = ',num2str(cell2mat(Coils0(3,3)))]);
% disp(['Io2 = ',num2str(VCI2(4)),'; Io20 = ',num2str(cell2mat(Coils0(4,3)))]);
% disp(['Io1 = ',num2str(VCI2(5)),'; Io10 = ',num2str(cell2mat(Coils0(5,3)))]);

disp(['plasma current = ',num2str(tmpI1),'; Ip0 = ',num2str(IDLP0(3,end))]);
