function [RHS] = p2ga(Coils,GR,GZ)
% ---- 选择的R和Z范围需要包含L线圈
% ---------- p DCoil 2 g ----------

cdr1=CRD(1,:);cdr2=CRD(end,:);cdr3=CRD(:,1)';cdr4=CRD(:,end)';
cdz1=CZD(1,:);cdz2=CZD(end,:);cdz3=CZD(:,1)';cdz4=CZD(:,end)';
cdr=[cdr1,cdr2,cdr3,cdr4];
cdz=[cdz1,cdz2,cdz3,cdz4];
cdb=[cdr;cdz];

[R,Z] = ndgrid(rr,zz);
RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);
[RHSD,~] = inpoly2([RCol ZCol],cdb');
RHSD = reshape(RHSD,length(rr),length(zz));
RHSD = RHSD*CJ(2)/sum(RHSD==1,'all'); % 每个网格的电流大小

% ---------- p LCoil 2 g ---------
if ( max(zz) > max(max(CZL)) )
    clr1=CRL(1,:);clr2=CRL(end,:);clr3=CRL(:,1)';clr4=CRL(:,end)';
    clz1=CZL(1,:);clz2=CZL(end,:);clz3=CZL(:,1)';clz4=CZL(:,end)';
    clr=[clr1,clr2,clr3,clr4];
    clz=[clz1,clz2,clz3,clz4];
    clb=[clr;clz];
    
    [RHSL,~] = inpoly2([RCol ZCol],clb');
    RHSL = reshape(RHSL,length(rr),length(zz));
    RHSL = RHSL*CJ(1)/sum(RHSL==1,'all');%除ds是电流密度
    
    RHS = RHSD + RHSL;
    
else
    
    RHS = RHSD;
end
end

% --------------------- Orgional Code ----------------------------------- %
% % ---- 选择的R和Z范围需要包含L线圈
% % ---------- p DCoil 2 g ----------
% 
% cdr1=CRD(1,:);cdr2=CRD(end,:);cdr3=CRD(:,1)';cdr4=CRD(:,end)';
% cdz1=CZD(1,:);cdz2=CZD(end,:);cdz3=CZD(:,1)';cdz4=CZD(:,end)';
% cdr=[cdr1,cdr2,cdr3,cdr4];
% cdz=[cdz1,cdz2,cdz3,cdz4];
% cdb=[cdr;cdz];
% 
% [R,Z] = ndgrid(rr,zz);
% RCol=reshape(R,[],1);ZCol=reshape(Z,[],1);
% [RHSD,~] = inpoly2([RCol ZCol],cdb');
% RHSD = reshape(RHSD,length(rr),length(zz));
% RHSD = RHSD*CJ(2)/sum(RHSD==1,'all');
% 
% % ---------- p LCoil 2 g ---------
% 
% clr1=CRL(1,:);clr2=CRL(end,:);clr3=CRL(:,1)';clr4=CRL(:,end)';
% clz1=CZL(1,:);clz2=CZL(end,:);clz3=CZL(:,1)';clz4=CZL(:,end)';
% clr=[clr1,clr2,clr3,clr4];
% clz=[clz1,clz2,clz3,clz4];
% clb=[clr;clz];
% 
% [RHSL,~] = inpoly2([RCol ZCol],clb');
% RHSL = reshape(RHSL,length(rr),length(zz));
% RHSL = RHSL*CJ(1)/sum(RHSL==1,'all');
% 
% 
% RHS = RHSD + RHSL;



