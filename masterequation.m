function rhodot = masterequation(t,rho,~,dim,t1,t2,Hzero,Hv,omega_0,Omega,delta,alpha,Gamma_0r,Gamma_0a,...
    Gamma_1r,Gamma_1a,Gamma_ar,Gamma_aa,gamma_r,s0,s1,sr,sa,plus,minus,decay0r_1,decay1r_1,decayar_1,decay0a_1,decay1a_1,decayaa_1,...
    decay0r_2,decay1r_2,decayar_2,decay0a_2,decay1a_2,decayaa_2)

dm = reshape(rho,dim^2,dim^2);

detuning = (2*delta*sin((omega_0*t)/2)^2)/omega_0;

HM = Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*minus' + (Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*minus')';
HP = Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*plus' + (Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*plus')';
H0 = Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*s0' + (Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*s0')';
H1 = Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*s1' + (Omega/2*exp(1i*(alpha+pi/2))*exp(-1i*detuning)*sr*s1')';

URP_12 = kron(HM, eye(4)) + kron(eye(4), H1) + Hv;
URP_23 = kron(H1, eye(4)) + kron(eye(4), HM) + Hv;

Hami = URP_12.*(t>0 & t<=t1) + Hzero.*(t>t1 & t<=t1+t2)...
     + URP_23.*(t>t1+t2 & t<=2*t1+t2) + Hzero.*(t>2*t1+1*t2 & t<=2*t1+2*t2);

%%
gamma_0r = 0.*(t>0 & t<=t1) + Gamma_0r.*(t>t1 & t<=t1+t2)...
    + gamma_r/2.*(t>t1+t2 & t<=2*t1+t2) + Gamma_0r.*(t>2*t1+1*t2 & t<=2*t1+2*t2);
gamma_1r = 0.*(t>0 & t<=t1) + Gamma_1r.*(t>t1 & t<=t1+t2)...
    + gamma_r/2.*(t>t1+t2 & t<=2*t1+t2) + Gamma_1r.*(t>2*t1+1*t2 & t<=2*t1+2*t2);
gamma_ar = 0.*(t>0 & t<=t1) + Gamma_ar.*(t>t1 & t<=t1+t2)...
    + gamma_r/2.*(t>t1+t2 & t<=2*t1+t2) + Gamma_ar.*(t>2*t1+1*t2 & t<=2*t1+2*t2);
gamma_0a = 0.*(t>0 & t<=t1) + Gamma_0a.*(t>t1 & t<=t1+t2)...
    + 0.*(t>t1+t2 & t<=2*t1+t2) + Gamma_0a.*(t>2*t1+1*t2 & t<=2*t1+2*t2);
gamma_1a = 0.*(t>0 & t<=t1) + Gamma_1a.*(t>t1 & t<=t1+t2)...
    + 0.*(t>t1+t2 & t<=2*t1+t2) + Gamma_1a.*(t>2*t1+1*t2 & t<=2*t1+2*t2);
gamma_aa = 0.*(t>0 & t<=t1) + Gamma_aa.*(t>t1 & t<=t1+t2)...
    + 0.*(t>t1+t2 & t<=2*t1+t2) + Gamma_aa.*(t>2*t1+1*t2 & t<=2*t1+2*t2);

temp0r_1 = 2*decay0r_1*dm*decay0r_1'-decay0r_1'*decay0r_1*dm-dm*decay0r_1'*decay0r_1;
temp1r_1 = 2*decay1r_1*dm*decay1r_1'-decay1r_1'*decay1r_1*dm-dm*decay1r_1'*decay1r_1;
tempar_1 = 2*decayar_1*dm*decayar_1'-decayar_1'*decayar_1*dm-dm*decayar_1'*decayar_1;
temp0a_1 = 2*decay0a_1*dm*decay0a_1'-decay0a_1'*decay0a_1*dm-dm*decay0a_1'*decay0a_1;
temp1a_1 = 2*decay1a_1*dm*decay1a_1'-decay1a_1'*decay1a_1*dm-dm*decay1a_1'*decay1a_1;
tempaa_1 = 2*decayaa_1*dm*decayaa_1'-decayaa_1'*decayaa_1*dm-dm*decayaa_1'*decayaa_1;

temp0r_2 = 2*decay0r_2*dm*decay0r_2'-decay0r_2'*decay0r_2*dm-dm*decay0r_2'*decay0r_2;
temp1r_2 = 2*decay1r_2*dm*decay1r_2'-decay1r_2'*decay1r_2*dm-dm*decay1r_2'*decay1r_2;
tempar_2 = 2*decayar_2*dm*decayar_2'-decayar_2'*decayar_2*dm-dm*decayar_2'*decayar_2;
temp0a_2 = 2*decay0a_2*dm*decay0a_2'-decay0a_2'*decay0a_2*dm-dm*decay0a_2'*decay0a_2;
temp1a_2 = 2*decay1a_2*dm*decay1a_2'-decay1a_2'*decay1a_2*dm-dm*decay1a_2'*decay1a_2;
tempaa_2 = 2*decayaa_2*dm*decayaa_2'-decayaa_2'*decayaa_2*dm-dm*decayaa_2'*decayaa_2;

out = -1i*Hami*dm+1i*dm*Hami...
    + 1/2*gamma_0r*temp0r_1 + 1/2*gamma_1r*temp1r_1 + 1/2*gamma_ar*tempar_1...
    + 1/2*gamma_0a*temp0a_1 + 1/2*gamma_1a*temp1a_1 + 1/2*gamma_aa*tempaa_1...
    + 1/2*gamma_0r*temp0r_2 + 1/2*gamma_1r*temp1r_2 + 1/2*gamma_ar*tempar_2...
    + 1/2*gamma_0a*temp0a_2 + 1/2*gamma_1a*temp1a_2 + 1/2*gamma_aa*tempaa_2;

rhodot = reshape(out,(dim^2)^2,1);