function [G, Iion, Ivec, dX, Iall] = model_ORd11_cable_alex_ORd2021(t,X,p,S,gendertype, par_SA)

% Units:
%
% time in msec
% voltage in mV
% current in uA
% concentration in mM
% conductance in mS
% capacitance in uF (uF*mV/ms = uA)

if ~isfield(p,'dt')
    p.dt = 1;
end
dt = p.dt;
Nm = p.Npatches;

INa_Multiplier = par_SA(1)*ones(Nm,1);   %'GNa' 
INaL_Multiplier = par_SA(2)*ones(Nm,1);  %'GNaL'
Ito_Multiplier = par_SA(3)*ones(Nm,1);   %'Gto'
ICaL_Multiplier = par_SA(4)*ones(Nm,1);  %'PCa'
IKr_Multiplier = par_SA(5)*ones(Nm,1);   %'GKr'
IKs_Multiplier = par_SA(6)*ones(Nm,1);   %'GKs'
IK1_Multiplier = par_SA(7)*ones(Nm,1);   %'GK1'
INaCa_Multiplier = par_SA(8)*ones(Nm,1); %'Gncx'
INaK_Multiplier = par_SA(9)*ones(Nm,1);  %'Pnak'
IKB_Multiplier = par_SA(10)*ones(Nm,1);  %'GKb'
INaB_Multiplier = par_SA(11)*ones(Nm,1); %'PNab'
ICaB_Multiplier = par_SA(12)*ones(Nm,1); %'PCab'
IpCa_Multiplier = par_SA(13)*ones(Nm,1); %'GpCa'
Jup_Multiplier = par_SA(14)*ones(Nm,1);  %'fSERCA'
Jrel_Multiplier = par_SA(15)*ones(Nm,1); %'fRyR'
Jleak_Multiplier = par_SA(16)*ones(Nm,1);%'fIleak'

v = X(1:Nm);  % voltages
G = X(Nm+1:end); % gating variables and ionic concentrations
dX = zeros(length(v)+length(G),1);

%extracellular ionic concentrations, mM

nao = S(1:Nm);   % extracellular Na
ko = S(Nm+1:2*Nm); % extracellular Ko
cao = S(2*Nm+1:3*Nm); % extracellular Ca

% additional scaling factors
fSERCA = p.fSERCA;
fRyR = p.fRyR;
ftauhL = p.ftauhL;
fCaMKa = p.fCaMKa;
fIleak = p.fIleak;
fJrel = p.fJrel;

celltype = p.celltype; %endo = 0, epi = 1, M = 2

%physical constants
R=8314.0;
T=310.0;
F=96485.0;  % C/mol

%cell geometry -
L=p.L/1e4;  % cm
rad=p.r/1e4;  % cm
vcell=1000*pi*rad*rad*L;  % uL
Ageo=2*pi*rad*rad+2*pi*rad*L;  % cm^2
Rcg = 2;
Acap=Rcg*Ageo;  % capacitive area, uF
vmyo=0.68*vcell;
vnsr=0.0552*vcell;
vjsr=0.0048*vcell;
vss=0.02*vcell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%give names to the state vector values
% v=X(1);
nai=X(1*Nm+1:2*Nm);
nass=X(2*Nm+1:3*Nm);
ki=X(3*Nm+1:4*Nm);
kss=X(4*Nm+1:5*Nm);
cai=X(5*Nm+1:6*Nm);
cass=X(6*Nm+1:7*Nm);
cansr=X(7*Nm+1:8*Nm);
cajsr=X(8*Nm+1:9*Nm);
m=X(9*Nm+1:10*Nm);
hf=X(10*Nm+1:11*Nm);
hs=X(11*Nm+1:12*Nm);
j=X(12*Nm+1:13*Nm);
hsp=X(13*Nm+1:14*Nm);
jp=X(14*Nm+1:15*Nm);
mL=X(15*Nm+1:16*Nm);
hL=X(16*Nm+1:17*Nm);
hLp=X(17*Nm+1:18*Nm);
a=X(18*Nm+1:19*Nm);
iF=X(19*Nm+1:20*Nm);
iS=X(20*Nm+1:21*Nm);
ap=X(21*Nm+1:22*Nm);
iFp=X(22*Nm+1:23*Nm);
iSp=X(23*Nm+1:24*Nm);
d=X(24*Nm+1:25*Nm);
ff=X(25*Nm+1:26*Nm);
fs=X(26*Nm+1:27*Nm);
fcaf=X(27*Nm+1:28*Nm);
fcas=X(28*Nm+1:29*Nm);
jca=X(29*Nm+1:30*Nm);
nca=X(30*Nm+1:31*Nm);
ffp=X(31*Nm+1:32*Nm);
fcafp=X(32*Nm+1:33*Nm);
xrf=X(33*Nm+1:34*Nm);
xrs=X(34*Nm+1:35*Nm);
xs1=X(35*Nm+1:36*Nm);
xs2=X(36*Nm+1:37*Nm);
xk1=X(37*Nm+1:38*Nm);
Jrelnp=X(38*Nm+1:39*Nm);
Jrelp=X(39*Nm+1:40*Nm);
CaMKt=X(40*Nm+1:41*Nm);

% 41 state variables, including Vm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CaMK constants
KmCaMK=0.15;

aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;
KmCaM=0.0015;
%update CaMK
CaMKb=CaMKo*(1.0-CaMKt)./(1.0+KmCaM./cass);
CaMKa=fCaMKa*(CaMKb+CaMKt);
dCaMKt=aCaMK*CaMKb.*(CaMKb+CaMKt)-bCaMK*CaMKt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reversal potentials, mV
ENa=(R*T/F)*log(nao./nai);
EK=(R*T/F)*log(ko./ki);
PKNa=0.01833;
EKs=(R*T/F)*log((ko+PKNa*nao)./(ki+PKNa*nai));

%convenient shorthand calculations
vffrt=v*F*F/(R*T);
vfrt=v*F/(R*T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currents calculated in uA/uF
%calculate INa
mss=1.0./(1.0+exp((-(v+39.57))/9.871));
tm=1.0./(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
% dm=(mss-m)./tm;
hss=1.0./(1+exp((v+82.90)/6.086));
thf=1.0./(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
ths=1.0./(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
Ahf=0.99;
Ahs=1.0-Ahf;
% dhf=(hss-hf)./thf;
% dhs=(hss-hs)./ths;
h=Ahf*hf+Ahs*hs;
jss=hss;
tj=2.038+1.0./(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
% dj=(jss-j)./tj;
hssp=1.0./(1+exp((v+89.1)/6.086));
thsp=3.0*ths;
% dhsp=(hssp-hsp)./thsp;
hp=Ahf*hf+Ahs*hsp;
tjp=1.46*tj;
% djp=(jss-jp)./tjp;
GNa=75;
fINap=(1.0./(1.0+KmCaMK./CaMKa));
INa=p.f_I(p.iina,:)'.*INa_Multiplier.*GNa.*(v-ENa).*m.^3.*((1.0-fINap).*h.*j+fINap.*hp.*jp);

%calculate INaL
mLss=1.0./(1.0+exp((-(v+42.85))/5.264));
tmL=tm;
% dmL=(mLss-mL)./tmL;
hLss=1.0./(1.0+exp((v+87.61)/7.488));
thL=ftauhL*200.0;
% dhL=(hLss-hL)./thL;
hLssp=1.0./(1.0+exp((v+93.81)/7.488));
thLp=3.0*thL;
% dhLp=(hLssp-hLp)./thLp;
GNaL=0.0075*ones(Nm,1);
GNaL(celltype==1) = 0.0075*0.6;
% if celltype(1)==1
%     GNaL=GNaL*0.6;
% end
fINaLp=(1.0./(1.0+KmCaMK./CaMKa));
INaL=p.f_I(p.iinal,:)'.*INaL_Multiplier.*GNaL.*(v-ENa).*mL.*((1.0-fINaLp).*hL+fINaLp.*hLp);

%calculate Ito
ass=1.0./(1.0+exp((-(v-14.34))/14.82));
ta=1.0515./(1.0./(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5./(1.0+exp((v+100.0)/29.3814)));
% da=(ass-a)./ta;
iss=1.0./(1.0+exp((v+43.94)/5.711));
delta_epi = ones(Nm,1);
tmp = (celltype ==1);
delta_epi(tmp) = 1.0-(0.95./(1.0+exp((v(tmp)+70.0)/5.0)));
% if celltype(1)==1
%     delta_epi=1.0-(0.95./(1.0+exp((v+70.0)/5.0)));
% else
%     delta_epi=1.0;
% end
tiF=4.562+1./(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
tiS=23.62+1./(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
tiF=tiF.*delta_epi;
tiS=tiS.*delta_epi;
AiF=1.0./(1.0+exp((v-213.6)/151.2));
AiS=1.0-AiF;
% diF=(iss-iF)./tiF;
% diS=(iss-iS)./tiS;

scaleItos = ones(Nm,1);
if gendertype == 1
    scaleItos(celltype==1) = 0.6;
    scaleItos(celltype==0) = 1;
elseif gendertype == 2
    scaleItos(celltype==1) = 0.26;
    scaleItos(celltype==0) = 0.64;
else 
    scaleItos(celltype==1) = 1;
    scaleItos(celltype==0) = 1;

end

i=AiF.*iF+AiS.*iS.*scaleItos;
assp=1.0./(1.0+exp((-(v-24.34))/14.82));
% dap=(assp-ap)./ta;
dti_develop=1.354+1.0e-4./(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
dti_recover=1.0-0.5./(1.0+exp((v+70.0)/20.0));
tiFp=dti_develop.*dti_recover.*tiF;
tiSp=dti_develop.*dti_recover.*tiS;
% diFp=(iss-iFp)./tiFp;
% diSp=(iss-iSp)./tiSp;

ip=AiF.*iFp+AiS.*iSp.*scaleItos;
Gto=0.02*ones(Nm,1);
Gto(celltype ==1) = 0.02*2;  %Scale by 2 instead of 4
%Gto(celltype(1)==2) = 0.02*4;
% if celltype(1)==1
%     Gto=Gto*2.0;
% elseif celltype(1)==2
%     Gto=Gto*4.0;
% end
fItop=(1.0./(1.0+KmCaMK./CaMKa));
Ito=p.f_I(p.iito,:)'.*Ito_Multiplier.*Gto.*(v-EK).*((1.0-fItop).*a.*i+fItop.*ap.*ip);

%calculate ICaL, ICaNa, ICaK
dss=1.0./(1.0+exp((-(v+3.940))/4.230));
td=0.6+1.0./(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
% dd=(dss-d)./td;
fss=1.0./(1.0+exp((v+19.58)/3.696));
tff=7.0+1.0./(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
tfs=1000.0+1.0./(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
Aff=0.6;
Afs=1.0-Aff;
% dff=(fss-ff)./tff;
% dfs=(fss-fs)./tfs;
f=Aff.*ff+Afs.*fs;
fcass=fss;
tfcaf=7.0+1.0./(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
tfcas=100.0+1.0./(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
Afcaf=0.3+0.6./(1.0+exp((v-10.0)/10.0));
Afcas=1.0-Afcaf;
% dfcaf=(fcass-fcaf)./tfcaf;
% dfcas=(fcass-fcas)./tfcas;
fca=Afcaf.*fcaf+Afcas.*fcas;
tjca=75.0;
% djca=(fcass-jca)./tjca;
tffp=2.5*tff;
% dffp=(fss-ffp)./tffp;
fp=Aff.*ffp+Afs.*fs;
tfcafp=2.5*tfcaf;
% dfcafp=(fcass-fcafp)./tfcafp;
fcap=Afcaf.*fcafp+Afcas.*fcas;
Kmn=0.002;
k2n=1000.0;
km2n=jca*1.0;
anca=1.0./(k2n./km2n+(1.0+Kmn./cass).^4.0);
dnca=anca.*k2n-nca.*km2n;
PhiCaL=4.0*vffrt.*(cass.*exp(2.0*vfrt)-0.341*cao)./(exp(2.0*vfrt)-1.0);
PhiCaNa=1.0*vffrt.*(0.75*nass.*exp(1.0*vfrt)-0.75*nao)./(exp(1.0*vfrt)-1.0);
PhiCaK=1.0*vffrt.*(0.75*kss.*exp(1.0*vfrt)-0.75*ko)./(exp(1.0*vfrt)-1.0);
zca=2.0;
PCa=0.0001*ones(Nm,1).*ICaL_Multiplier;
PCa(celltype==1) = 0.0001*1.2;
%PCa(celltype ==2) = 0.0001*2.5;
% if celltype(1)==1
%     PCa=PCa*1.2;
% elseif celltype(1)==2
%     PCa=PCa*2.5;
% end
PCap=1.1*PCa;
PCaNa=0.00125*PCa;
PCaK=3.574e-4*PCa;
PCaNap=0.00125*PCap;
PCaKp=3.574e-4*PCap;
fICaLp=(1.0./(1.0+KmCaMK./CaMKa));
ICaL=p.f_I(p.iical,:)'.*((1.0-fICaLp).*PCa.*PhiCaL.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCap.*PhiCaL.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca));
ICaNa=p.f_I(p.iical,:)'.*((1.0-fICaLp).*PCaNa.*PhiCaNa.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCaNap.*PhiCaNa.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca));
ICaK=p.f_I(p.iical,:)'.*((1.0-fICaLp).*PCaK.*PhiCaK.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp.*PCaKp.*PhiCaK.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca));

%calculate IKr
xrss=1.0./(1.0+exp((-(v+8.337))/6.789));
txrf=12.98+1.0./(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
txrs=1.865+1.0./(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
Axrf=1.0./(1.0+exp((v+54.81)/38.21));
Axrs=1.0-Axrf;
% dxrf=(xrss-xrf)./txrf;
% dxrs=(xrss-xrs)./txrs;
xr=Axrf.*xrf+Axrs.*xrs;
rkr=1.0./(1.0+exp((v+55.0)/75.0))*1.0./(1.0+exp((v-10.0)/30.0));
% GKr=0.046*ones(Nm,1);
% GKr(celltype(1)==1) = 0.046*1.3;
% GKr(celltype(1)==2) = 0.046*0.8;
% if celltype(1)==1
%     GKr=GKr*1.3;
% elseif celltype(1)==2
%     GKr=GKr*0.8;
% end

%AP1_gKr = 0.04; AP165_gKr = 0.06; %1.5b

if gendertype == 1
    GKr= p.GKr_male;%linspace(AP1_gKr, AP165_gKr, Nm )';
    %GKr= linspace(0.02, 0.05, Nm )';
    %GKr= linspace(0.046, 0.05, Nm )';
    % GKr=0.046*ones(Nm,1);
    % GKr(celltype==1) = GKr(celltype==1)*1.09;
    % GKr(celltype==0) = GKr(celltype==0)*1;
elseif gendertype == 2
    GKr= p.GKr_female;%linspace(0.79*AP1_gKr, 0.875*AP165_gKr, Nm )';
    %GKr= linspace(0.02*0.875, 0.02*0.875*2.52, Nm )';
    %GKr= linspace(0.036, 0.042, Nm)';
    % GKr=0.046*ones(Nm,1);
    % GKr(celltype==1) = GKr(celltype==1)*0.875;
    % GKr(celltype==0) = GKr(celltype==0)*0.79;
else 
    GKr=0.046*ones(Nm,1);
    GKr(celltype==1) = 0.046*1.3;  %0.06
    GKr(celltype==2) = 0.046*0.8;  %0.36
end


IKr=p.f_I(p.iikr,:)'.*IKr_Multiplier.*GKr.*sqrt(ko/5.4).*xr.*rkr.*(v-EK);

%calculate IKs
% if celltype(1)==1
%     GKs=GKs*1.4;
% end
GKs=0.0034*ones(Nm,1);
xs1ss=1.0./(1.0+exp((-(v+11.60))/8.932));
txs1=817.3+1.0./(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
xs2ss=xs1ss;
txs2=1.0./(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));

if gendertype == 1
    xs1ss(celltype==1)=1.0./(1.0+exp((-(v(celltype==1)+11.60))/8.932)*1.04);
    txs1(celltype==1)=817.3+1.0./(2.326e-4*exp((v(celltype==1)+48.28)/17.80)*1.04+0.001292*exp((-(v(celltype==1)+210.0))/230.0)*1.04);
    xs2ss=xs1ss;
    txs2(celltype==1)=1.0./(0.01*exp((v(celltype==1)-50.0)/20.0)*1.04+0.0193*exp((-(v(celltype==1)+66.54))/31.0)*1.04);
    GKs(celltype==1) = 0.0034*1.04;
elseif gendertype == 2
    xs1ss(celltype==1)=1.0./(1.0+exp((-(v(celltype==1)+11.60))/8.932)*0.87);
    txs1(celltype==1)=817.3+1.0./(2.326e-4*exp((v(celltype==1)+48.28)/17.80)*0.87+0.001292*exp((-(v(celltype==1)+210.0))/230.0)*0.87);
    txs2(celltype==1)=1.0./(0.01*exp((v(celltype==1)-50.0)/20.0)*0.87+0.0193*exp((-(v(celltype==1)+66.54))/31.0)*0.87);
    GKs(celltype==1) = 0.0034*0.87;

    xs1ss(celltype==0)=1.0./(1.0+exp((-(v(celltype==0)+11.60))/8.932)*0.83);
    txs1(celltype==0)=817.3+1.0./(2.326e-4*exp((v(celltype==0)+48.28)/17.80)*0.83+0.001292*exp((-(v(celltype==0)+210.0))/230.0)*0.83);
    xs2ss=xs1ss;
    txs2(celltype==0)=1.0./(0.01*exp((v(celltype==0)-50.0)/20.0)*0.83+0.0193*exp((-(v(celltype==0)+66.54))/31.0)*0.83);
    GKs(celltype==0) = 0.0034*0.83;
else
    GKs(celltype ==1) = 0.0034*1.4;
end

KsCa=1.0+0.6./(1.0+(3.8e-5./cai).^1.4);
IKs=p.f_I(p.iiks,:)'.*IKs_Multiplier.*GKs.*KsCa.*xs1.*xs2.*(v-EKs);

% calculate IK1
xk1ss=1.0./(1.0+exp(-(v+2.5538*ko+144.59)./(1.5692*ko+3.8115)));
txk1=122.2./(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
% dxk1=(xk1ss-xk1)./txk1;
rk1=1.0./(1.0+exp((v+105.8-2.6*ko)./9.493));
GK1=0.1908*ones(Nm,1);

% if celltype(1)==1
%     GK1=GK1*1.2;
% elseif celltype(1)==2
%     GK1=GK1*1.3;
% end
if gendertype == 1
    GK1(celltype==1) = 0.1908*0.98;
    GK1(celltype==0) = 0.1908*1;
elseif gendertype == 2
    GK1(celltype==1) = 0.1908*0.74;
    GK1(celltype==0) = 0.1908*0.86;
else
    GK1(celltype==1) = 0.1908*1.2;
    GK1(celltype==2) = 0.1908*1.3;
end

IK1=p.f_I(p.iik1,:)'.*IK1_Multiplier.*GK1.*sqrt(ko).*rk1.*xk1.*(v-EK);

%calculate INaCa_i
kna1=15.0;
kna2=5.0;
kna3=88.12;
kasymm=12.5;
wna=6.0e4;
wca=6.0e4;
wnaca=5.0e3;
kcaon=1.5e6;
kcaoff=5.0e3;
qna=0.5224;
qca=0.1670;
hca=exp((qca*v*F)/(R*T));
hna=exp((qna*v*F)/(R*T));
h1=1+nai./kna3.*(1+hna);
h2=(nai.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+nai./kna1.*(1+nai./kna2);
h5=nai.*nai./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+nao./kna3.*(1.0+1.0./hna);
h8=nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+nao./kna1.*(1.0+nao./kna2);
h11=nao.*nao./(h10.*kna1.*kna2);
h12=1.0./h10;
k1=h12.*cao.*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*cai.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./cai).^2.0);
zna=1.0;
JncxNa=3.0*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
Gncx=0.0008*ones(Nm,1);
Gncx(celltype ==1) = 0.0008*1.1;  
% Gncx(celltype(1)==2) = 0.0008*1.4;
% if celltype(1)==1
%     Gncx=Gncx*1.1;   
% elseif celltype(1)==2
%     Gncx=Gncx*1.4;
%end
if gendertype==1
    scaleNaCas=1;
elseif gendertype==2
    scaleNaCas=1.15;
else
    scaleNaCas=1;
end


INaCa_i=p.f_I(p.iinaca_i,:)'.*0.8.*INaCa_Multiplier.*Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa)*scaleNaCas;

%calculate INaCa_ss
h1=1+nass./kna3.*(1+hna);
h2=(nass.*hna)./(kna3.*h1);
h3=1.0./h1;
h4=1.0+nass./kna1.*(1+nass./kna2);
h5=nass.*nass./(h4.*kna1.*kna2);
h6=1.0./h4;
h7=1.0+nao./kna3.*(1.0+1.0./hna);
h8=nao./(kna3.*hna.*h7);
h9=1.0./h7;
h10=kasymm+1.0+nao./kna1.*(1+nao./kna2);
h11=nao.*nao./(h10.*kna1.*kna2);
h12=1.0./h10;
k1=h12.*cao.*kcaon;
k2=kcaoff;
k3p=h9.*wca;
k3pp=h8.*wnaca;
k3=k3p+k3pp;
k4p=h3.*wca./hca;
k4pp=h2.*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6.*cass.*kcaon;
k7=h5.*h2.*wna;
k8=h8.*h11.*wna;
x1=k2.*k4.*(k7+k6)+k5.*k7.*(k2+k3);
x2=k1.*k7.*(k4+k5)+k4.*k6.*(k1+k8);
x3=k1.*k3.*(k7+k6)+k8.*k6.*(k2+k3);
x4=k2.*k8.*(k4+k5)+k3.*k5.*(k1+k8);
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0./(1.0+(KmCaAct./cass).^2.0);
JncxNa=3.0*(E4.*k7-E1.*k8)+E3.*k4pp-E2.*k3pp;
JncxCa=E2.*k2-E1.*k1;
INaCa_ss=p.f_I(p.iinaca_ss,:)'.*0.2.*INaCa_Multiplier.*Gncx.*allo.*(zna.*JncxNa+zca.*JncxCa).*scaleNaCas;

%calculate INaK
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;
Knai0=9.073;
Knao0=27.78;
delta=-0.1550;
Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
Kki=0.5;
Kko=0.3582;
MgADP=0.05;
MgATP=9.8;
Kmgatp=1.698e-7;
H=1.0e-7;
eP=4.2;
Khp=1.698e-7;
Knap=224.0;
Kxkur=292.0;
P=eP./(1.0+H/Khp+nai./Knap+ki./Kxkur);
a1=(k1p*(nai./Knai).^3.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
b1=k1m*MgADP;
a2=k2p;
b2=(k2m*(nao./Knao).^3.0)./((1.0+nao./Knao).^3.0+(1.0+ko./Kko).^2.0-1.0);
a3=(k3p*(ko/Kko).^2.0)./((1.0+nao./Knao).^3.0+(1.0+ko./Kko).^2.0-1.0);
b3=(k3m*P.*H)./(1.0+MgATP./Kmgatp);
a4=(k4p*MgATP/Kmgatp)./(1.0+MgATP./Kmgatp);
b4=(k4m.*(ki./Kki).^2.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
x1=a4.*a1.*a2+b2.*b4.*b3+a2.*b4.*b3+b3.*a1.*a2;
x2=b2.*b1.*b4+a1.*a2.*a3+a3.*b1.*b4+a2.*a3.*b4;
x3=a2.*a3.*a4+b3.*b2.*b1+b2.*b1.*a4+a3.*a4.*b1;
x4=b4.*b3.*b2+a3.*a4.*a1+b2.*a4.*a1+b3.*b2.*a1;
E1=x1./(x1+x2+x3+x4);
E2=x2./(x1+x2+x3+x4);
E3=x3./(x1+x2+x3+x4);
E4=x4./(x1+x2+x3+x4);
zk=1.0;
JnakNa=3.0*(E1.*a3-E2.*b3);
JnakK=2.0*(E4.*b1-E3.*a1);
Pnak=30*ones(Nm,1);

% if celltype(1)==1
%     Pnak=Pnak*0.9;
% elseif celltype(1)==2
%     Pnak=Pnak*0.7;
% end
%if gendertype == 1
    Pnak(celltype==1) = 30*0.94;
    Pnak(celltype==0) = 30*1;
%elseif gendertype == 2
    %Pnak(celltype==1) = 30*0.7;
    %Pnak(celltype==0) = 30*0.79;
%else
    %Pnak(celltype(1)==1) = 30*0.9;
    %Pnak(celltype(1)==2) = 30*0.7;
%end

INaK=p.f_I(p.iinak,:)'.*INaK_Multiplier.*Pnak.*(zna.*JnakNa+zk.*JnakK);

%calculate IKb
xkb=1.0./(1.0+exp(-(v-14.48)/18.34));
GKb=0.003*ones(Nm,1);
GKb(celltype ==1) = 0.003*0.6;
% if celltype(1)==1
%     GKb=GKb*0.6;
% end
IKb=p.f_I(p.iikb,:)'.*IKB_Multiplier.*GKb.*xkb.*(v-EK);

%calculate INab
PNab=3.75e-10;
INab=p.f_I(p.iinab,:)'.*INaB_Multiplier.*PNab.*vffrt.*(nai.*exp(vfrt)-nao)./(exp(vfrt)-1.0);

%calculate ICab
PCab=2.5e-8;
ICab=p.f_I(p.iicab,:)'.*ICaB_Multiplier.*PCab.*4.0.*vffrt.*(cai.*exp(2.0*vfrt)-0.341*cao)./(exp(2.0*vfrt)-1.0);

%calculate IpCa
GpCa=0.0005*ones(Nm,1);

if gendertype == 1
    GpCa(celltype==1) = 0.0005*0.88;
    GpCa(celltype==0) = 0.0005*1;
elseif gendertype == 2
    GpCa(celltype==1) = 0.0005*1.6;
    GpCa(celltype==0) = 0.0005*1.6;
end

IpCa=p.f_I(p.iipca,:)'.*IpCa_Multiplier.*GpCa.*cai./(0.0005+cai);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the stimulus current, Istim
% amp=80.0;
% duration=0.5;
% if t<=duration
%     Istim=amp;
% else
%     Istim=0.0;
% end

Istim = p.stim_amp*(mod(t,p.bcl)<p.stim_dur).*p.indstim;  % uA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% currents in uA/uF
% Iion, uA - need to scale by whole cell capacitance
Iion= p.Ctot*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab-Istim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate diffusion fluxes
JdiffNa=(nass-nai)/2.0;  % mM/msec
JdiffK=(kss-ki)/2.0;
Jdiff=(cass-cai)/0.2;

%calculate ryanodione receptor calcium induced calcium release from the jsr
bt=4.75;
a_rel=0.5*bt;
% ICaL_norm = ICaL./p.f_I(p.iical,:)';

Jrel_inf=fJrel*a_rel*(-ICaL)./(1.0+(1.5./cajsr).^8.0);
% tmp = (celltype(1)==2);
% Jrel_inf(tmp) = Jrel_inf(tmp)*1.7;
% if celltype(1)==2
%     Jrel_inf=Jrel_inf*1.7;
% end
tau_rel=bt./(1.0+0.0123./cajsr);

tau_rel(tau_rel<0.001) = 0.001;
% if tau_rel<0.001
%    tau_rel=0.001;
% end

% dJrelnp=(Jrel_inf-Jrelnp)./tau_rel;
btp=1.25*bt;
a_relp=0.5*btp;
Jrel_infp=a_relp*(-ICaL)./(1.0+(1.5./cajsr).^8.0);
%Jrel_infp(tmp) = Jrel_infp(tmp)*1.7;
% if celltype(1)==2
%     Jrel_infp=Jrel_infp*1.7;
% end
tau_relp=btp./(1.0+0.0123./cajsr);

tau_relp(tau_relp<0.001) = 0.001;

% if tau_relp<0.001
%    tau_relp=0.001;
% end

% dJrelp=(Jrel_infp-Jrelp)./tau_relp;
fJrelp=(1.0./(1.0+KmCaMK./CaMKa));
Jrel=Jrel_Multiplier.*fRyR.*((1.0-fJrelp).*Jrelnp+fJrelp.*Jrelp);

%calculate serca pump, ca uptake flux
Jupnp=0.004375*cai./(cai+0.00092);
Jupp=2.75*0.004375*cai./(cai+0.00092-0.00017);

% if celltype(1)==1
%     Jupnp=Jupnp*1.3;
%     Jupp=Jupp*1.3;
% end

%if gendertype == 1
    Jupnp(celltype==1) = Jupnp(celltype==1)*1.42;
    Jupp(celltype==1) = Jupp(celltype==1)*1.42;
    Jupnp(celltype==0) = Jupnp(celltype==0)*1;
    Jupp(celltype==0) = Jupp(celltype==0)*1;
%elseif gendertype == 2
    % Jupnp(celltype==1) = Jupnp(celltype==1)*1.97;
    % Jupp(celltype==1) = Jupp(celltype==1)*1.97;
    % Jupnp(celltype==0) = Jupnp(celltype==0)*1.15;
    % Jupp(celltype==0) = Jupp(celltype==0)*1.15;

% else
%     tmp = (celltype==1);
%     Jupnp(tmp) = Jupnp(tmp)*1.3;
%     Jupp(tmp) = Jupp(tmp)*1.3;
% 
% end

fJupp=(1.0./(1.0+KmCaMK./CaMKa));
Jleak=Jleak_Multiplier.*fIleak*0.0039375.*cansr/15.0;
Jup=Jup_Multiplier.*fSERCA.*((1.0-fJupp).*Jupnp+fJupp.*Jupp-Jleak);

%calculate tranlocation flux
Jtr=(cansr-cajsr)/100.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calcium buffer constants
cmdnmax=0.05*ones(Nm,1);

% if celltype(1)==1
%     cmdnmax=cmdnmax*1.3;
% end

if gendertype == 1
    cmdnmax(celltype==1) = 0.05*1.07;
    cmdnmax(celltype==0) = 0.05*1;
elseif gendertype == 2
    cmdnmax(celltype==1) = 0.05*1.41;
    cmdnmax(celltype==0) = 0.05*1.21;
else
    cmdnmax(celltype==1) = 0.05*1.3;
end

kmcmdn=0.00238;
trpnmax=0.07;
kmtrpn=0.0005;
BSRmax=0.047;
KmBSR=0.00087;
BSLmax=1.124;
KmBSL=0.0087;
csqnmax=10.0;
kmcsqn=0.8;

Ivec = nan(Nm*3,1);    % ionic currents, uA
Ivec(1:Nm) = p.Ctot*(INa+INaL+3*INaCa_i+3*INaK+INab + 3*INaCa_ss+ICaNa);      % Na+
Ivec(Nm+1:2*Nm) = p.Ctot*(Ito+IKr+IKs+IK1+IKb-Istim-2*INaK + ICaK); % K+
Ivec(2*Nm+1:3*Nm) = p.Ctot*(IpCa+ICab-2*INaCa_i + ICaL-2*INaCa_ss);      % Ca2+

% all ionic currents
[Ncurrents,~] = size(p.f_I);
Iall = nan(Nm*Ncurrents,1);

Iall(p.iina:Ncurrents:end) = p.Ctot*INa;
Iall(p.iinal:Ncurrents:end) = p.Ctot*INaL;
Iall(p.iito:Ncurrents:end) = p.Ctot*Ito;
Iall(p.iical:Ncurrents:end) = p.Ctot*ICaL;
Iall(p.iikr:Ncurrents:end) = p.Ctot*IKr;
Iall(p.iiks:Ncurrents:end) = p.Ctot*IKs;
Iall(p.iik1:Ncurrents:end) = p.Ctot*IK1;
Iall(p.iinaca_i:Ncurrents:end) = p.Ctot*INaCa_i;
Iall(p.iinaca_ss:Ncurrents:end) = p.Ctot*INaCa_ss;
Iall(p.iinak:Ncurrents:end) = p.Ctot*INaK;
Iall(p.iikb:Ncurrents:end) = p.Ctot*IKb;
Iall(p.iinab:Ncurrents:end) = p.Ctot*INab;
Iall(p.iicab:Ncurrents:end) = p.Ctot*ICab;
Iall(p.iipca:Ncurrents:end) = p.Ctot*IpCa;


% if any(imag(Ivec))
%     1;
% end
%update intracellular concentrations, using buffers for cai, cass, cajsr
%  units: (uA/uF)*uF / (C/mol*uL) = mM/msec
JintraNai = sum((nai(p.ind_tau_ip) - nai(p.ind_tau_im))./p.tau_Nai_mat,2);
JintraCai = sum((cai(p.ind_tau_ip) - cai(p.ind_tau_im))./p.tau_Cai_mat,2);
JintraKi = sum((ki(p.ind_tau_ip) - ki(p.ind_tau_im))./p.tau_Ki_mat,2);

dnai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo + JintraNai;
dnass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

dki=-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo + JintraKi;
dkss=-(ICaK)*Acap/(F*vss)-JdiffK;

Bcai=1.0./(1.0+cmdnmax.*kmcmdn./(kmcmdn+cai).^2.0+trpnmax.*kmtrpn./(kmtrpn+cai).^2.0);
dcai=Bcai.*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo + JintraCai);

Bcass=1.0./(1.0+BSRmax.*KmBSR./(KmBSR+cass).^2.0+BSLmax.*KmBSL./(KmBSL+cass).^2.0);
dcass=Bcass.*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

dcansr=Jup-Jtr*vjsr/vnsr;

Bcajsr=1.0./(1.0+csqnmax.*kmcsqn./(kmcsqn+cajsr).^2.0);
dcajsr=Bcajsr.*(Jtr-Jrel);

% output=[dv dnai dnass dki dkss dcai dcass dcansr dcajsr

% Forwad Euler, to update concentrations
% Rush-Larsen, to update gating variables
G(1:Nm) = nai + dt*dnai;
G(Nm+1:2*Nm) = nass + dt*dnass;
G(2*Nm+1:3*Nm) = ki + dt*dki;
G(3*Nm+1:4*Nm) = kss + dt*dkss;
G(4*Nm+1:5*Nm) = cai + dt*dcai;
G(5*Nm+1:6*Nm) = cass + dt*dcass;
G(6*Nm+1:7*Nm) = cansr + dt*dcansr;
G(7*Nm+1:8*Nm) = cajsr + dt*dcajsr;
G(8*Nm+1:9*Nm) =  mss - (mss - m).*exp(-dt./tm);
G(9*Nm+1:10*Nm) = hss - (hss - hf).*exp(-dt./thf);
G(10*Nm+1:11*Nm) = hss - (hss - hs).*exp(-dt./ths);
G(11*Nm+1:12*Nm) = jss - (jss - j).*exp(-dt./tj);
G(12*Nm+1:13*Nm) = hssp - (hssp - hsp).*exp(-dt./thsp);
G(13*Nm+1:14*Nm) = jss - (jss - jp).*exp(-dt./tjp);
G(14*Nm+1:15*Nm) = mLss - (mLss - mL).*exp(-dt./tmL);
G(15*Nm+1:16*Nm) = hLss - (hLss - hL).*exp(-dt./thL);
G(16*Nm+1:17*Nm) = hLssp - (hLssp - hLp).*exp(-dt./thLp);
G(17*Nm+1:18*Nm) = ass - (ass - a).*exp(-dt./ta);
G(18*Nm+1:19*Nm) = iss - (iss - iF).*exp(-dt./tiF);
G(19*Nm+1:20*Nm) = iss - (iss - iS).*exp(-dt./tiS);
G(20*Nm+1:21*Nm) = assp - (assp - ap).*exp(-dt./ta);
G(21*Nm+1:22*Nm) = iss - (iss - iFp).*exp(-dt./tiFp);
G(22*Nm+1:23*Nm) = iss - (iss - iSp).*exp(-dt./tiSp);
G(23*Nm+1:24*Nm) = dss - (dss - d).*exp(-dt./td);
G(24*Nm+1:25*Nm) = fss - (fss - ff).*exp(-dt./tff);
G(25*Nm+1:26*Nm) = fss - (fss - fs).*exp(-dt./tfs);
G(26*Nm+1:27*Nm) = fcass - (fcass - fcaf).*exp(-dt./tfcaf);
G(27*Nm+1:28*Nm) = fcass - (fcass - fcas).*exp(-dt./tfcas);
G(28*Nm+1:29*Nm) = fcass - (fcass - jca).*exp(-dt./tjca);
G(29*Nm+1:30*Nm) = nca + dt*dnca;
G(30*Nm+1:31*Nm) = fss - (fss - ffp).*exp(-dt./tffp);
G(31*Nm+1:32*Nm) = fcass - (fcass - fcafp).*exp(-dt./tfcafp);
G(32*Nm+1:33*Nm) = xrss - (xrss - xrf).*exp(-dt./txrf);
G(33*Nm+1:34*Nm) = xrss - (xrss - xrs).*exp(-dt./txrs);
G(34*Nm+1:35*Nm) = xs1ss - (xs1ss - xs1).*exp(-dt./txs1);
G(35*Nm+1:36*Nm) = xs2ss - (xs2ss - xs2).*exp(-dt./txs2);
G(36*Nm+1:37*Nm) = xk1ss - (xk1ss - xk1).*exp(-dt./txk1);
G(37*Nm+1:38*Nm) = Jrel_inf - (Jrel_inf - Jrelnp).*exp(-dt./tau_rel);
G(38*Nm+1:39*Nm) = Jrel_infp - (Jrel_infp - Jrelp).*exp(-dt./tau_relp);
G(39*Nm+1:40*Nm) = CaMKt + dt*dCaMKt;

if Nm == 1
    dX(1:Nm) = -Iion/p.Ctot;
    dX(Nm+1:2*Nm) = dnai;
    dX(2*Nm+1:3*Nm) = dnass;
    dX(3*Nm+1:4*Nm) = dki;
    dX(4*Nm+1:5*Nm) = dkss;
    dX(5*Nm+1:6*Nm) = dcai;
    dX(6*Nm+1:7*Nm) = dcass;
    dX(7*Nm+1:8*Nm) = dcansr;
    dX(8*Nm+1:9*Nm) = dcajsr;
    dX(9*Nm+1:10*Nm) =  (mss - m)./tm;
    dX(10*Nm+1:11*Nm) = (hss - hf)./thf;
    dX(11*Nm+1:12*Nm) = (hss - hs)./ths;
    dX(12*Nm+1:13*Nm) = (jss - j)./tj;
    dX(13*Nm+1:14*Nm) = (hssp - hsp)./thsp;
    dX(14*Nm+1:15*Nm) = (jss - jp)./tjp;
    dX(15*Nm+1:16*Nm) = (mLss - mL)./tmL;
    dX(16*Nm+1:17*Nm) = (hLss - hL)./thL;
    dX(17*Nm+1:18*Nm) = (hLssp - hLp)./thLp;
    dX(18*Nm+1:19*Nm) = (ass - a)./ta;
    dX(19*Nm+1:20*Nm) = (iss - iF)./tiF;
    dX(20*Nm+1:21*Nm) = (iss - iS)./tiS;
    dX(21*Nm+1:22*Nm) = (assp - ap)./ta;
    dX(22*Nm+1:23*Nm) =  (iss - iFp)./tiFp;
    dX(23*Nm+1:24*Nm) = (iss - iSp)./tiSp;
    dX(24*Nm+1:25*Nm) = (dss - d)./td;
    dX(25*Nm+1:26*Nm) = (fss - ff)./tff;
    dX(26*Nm+1:27*Nm) = (fss - fs)./tfs;
    dX(27*Nm+1:28*Nm) =  (fcass - fcaf)./tfcaf;
    dX(28*Nm+1:29*Nm) = (fcass - fcas)./tfcas;
    dX(29*Nm+1:30*Nm) = (fcass - jca)./tjca;
    dX(30*Nm+1:31*Nm) = dnca;
    dX(31*Nm+1:32*Nm) = (fss - ffp)./tffp;
    dX(32*Nm+1:33*Nm) = (fcass - fcafp)./tfcafp;
    dX(33*Nm+1:34*Nm) = (xrss - xrf)./txrf;
    dX(34*Nm+1:35*Nm) = (xrss - xrs)./txrs;
    dX(35*Nm+1:36*Nm) = (xs1ss - xs1)./txs1;
    dX(36*Nm+1:37*Nm) = (xs2ss - xs2)./txs2;
    dX(37*Nm+1:38*Nm) = (xk1ss - xk1)./txk1;
    dX(38*Nm+1:39*Nm) = (Jrel_inf - Jrelnp)./tau_rel;
    dX(39*Nm+1:40*Nm) = (Jrel_infp - Jrelp)./tau_relp;
    dX(40*Nm+1:41*Nm) = dCaMKt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


