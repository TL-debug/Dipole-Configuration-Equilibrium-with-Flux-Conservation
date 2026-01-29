function [SI,MDL] = SMI(VCRD,VCZD,R,Z,PSICD,PSICL,CJ)
% SMI 据公式计算D线圈自感 及 D和L线圈互感
% Self & Mult Introduce
% 1 D线圈自感
dr = R(2,1) - R(1,1);dz = Z(1,2) - Z(1,1);Ncd = size(VCRD,2)*size(VCZD,2);
[~,~,BBZ] = GradBB(PSICD,R,Z);
ZD = (max(VCZD) + min(VCZD))/2;
VRD = R(1:floor(VCRD(1)/dr) + 1,1)';
VZD = ZD*ones(size(VRD));
VBD = interp2(R',Z',BBZ',VRD,VZD);
SI = Ncd*2*pi*sum(VRD.*abs(VBD))*dr/CJ(2);


% 2 D和L线圈互感
% L线圈通电在D线圈的互感
[~,~,BBZ] = GradBB(PSICL,R,Z);
BBZ(1,:) = BBZ(2,:);
ZD = (max(VCZD) + min(VCZD))/2;
VRD = R(1:floor(VCRD(1)/dr) + 1,1)';
VZD = ZD*ones(size(VRD));
VBD = interp2(R',Z',BBZ',VRD,VZD);
MDL = Ncd*2*pi*sum(VRD.*abs(VBD))*dr/CJ(1);

% % L线圈平面的位置（差值中平面位置）
% [~,~,BBZ] = GradBB(PSICD,R,Z);
% BBZ(1,:) = BBZ(2,:);
% Ncl = size(VCRL,2)*size(VCZL,2);
% ZL = (max(VCZL) + min(VCZL))/2;
% VRL = R(1:floor(VCRL(1)/dr) + 1,1)';
% VZL = ZL*ones(size(VRL));
% VBL = interp2(R',Z',BBZ',VRL,VZL);
% MLD = Ncl*2*pi*sum(VRL.*abs(VBL))*dr/CJ(2);
end