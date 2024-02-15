clc;  clear all;  close all;
tic
dim = 4;
%% 设置算子
s0 = [1;0;0;0];   s1 = [0;1;0;0];   sr = [0;0;1;0];   sa = [0;0;0;1];
plus = (s0 + s1)/sqrt(2);    minus = (s0 - s1)/sqrt(2);
s0r = s0*sr';   s1r = s1*sr';   sar = sa*sr';
s0a = s0*sa';   s1a = s1*sa';   saa = sa*sa';
%% 设置decay算子
decay0r_1 = kron(s0r,eye(4)); decay1r_1 = kron(s1r,eye(4)); decayar_1 = kron(sar,eye(4));
decay0a_1 = kron(s0a,eye(4)); decay1a_1 = kron(s1a,eye(4)); decayaa_1 = kron(saa,eye(4));
decay0r_2 = kron(eye(4),s0r); decay1r_2 = kron(eye(4),s1r); decayar_2 = kron(eye(4),sar);
decay0a_2 = kron(eye(4),s0a); decay1a_2 = kron(eye(4),s1a); decayaa_2 = kron(eye(4),saa);
%% 定义稳定子
sigma_x = s0*s1' + s1*s0';   sigma_y = -1i*s0*s1' + 1i*s1*s0';   sigma_z = s0*s0' - s1*s1';
XZ = kron(sigma_x,sigma_z);
ZX = kron(sigma_z,sigma_x);
% 定义 Bell 态
Bell_00 = ( kron(s0,plus) + kron(s1,minus) )/sqrt(2);
Bell_01 = ( kron(s0,plus) - kron(s1,minus) )/sqrt(2);
Bell_10 = ( kron(s1,plus) + kron(s0,minus) )/sqrt(2);
Bell_11 = ( kron(s1,plus) - kron(s0,minus) )/sqrt(2);
%% 参数
lifetime_r = 155*10^-6;  lifetime_P1 = 26.2*10^-9;  lifetime_P2 = 27.7*10^-9;
gamma_r = 1/lifetime_r;    gamma_P1 = 1/lifetime_P1;    gamma_P2 = 1/lifetime_P2;
Omega_r = 0.2*gamma_P1;   Omega_a = 0.2*gamma_P2;
Gamma_0r = Omega_r^2/2/gamma_P1;  Gamma_1r = Omega_r^2/3/gamma_P1;  Gamma_ar = Omega_r^2/6/gamma_P1;
Gamma_0a = Omega_a^2/2/gamma_P2;  Gamma_1a = Omega_a^2/3/gamma_P2;  Gamma_aa = Omega_a^2/6/gamma_P2;

% Omega = 2*pi*1*10^6;   Urr = 2*pi*70*10^6;
% omega_0 = Urr/15; alpha = 49.506177;
Omega = 2*pi*1*10^6;   Urr = 2*pi*200*10^6;
omega_0 = -Urr/15; alpha = 49.506177;
delta = alpha*omega_0;
t1 = pi/(Omega*besselj(0,alpha));
t2 = 3*pi/(Omega_a^2/sqrt(gamma_P1^2+gamma_P2^2));
time = 2*t1+2*t2;
%% 初态 和 末态
psi_0 = kron(plus, plus);
% psi_0 = (Bell_00 + Bell_01 + Bell_10 + Bell_11)/sqrt(4);
psi_f = Bell_00;
dm1 = psi_0*psi_0';
fidm = psi_f*psi_f';

Hv = Urr*kron(sr*sr', sr*sr');

Hzero = zeros(dim^2,dim^2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1000;   timestep1 = time/N;   tlist = linspace(0,time,N);
options = odeset('reltol',1e-6,'abstol',1e-6);

for jj = 1:18
    [t,rho]=ode45('masterequation',tlist,dm1,options,dim,t1,t2,Hzero,Hv,omega_0,Omega,delta,alpha,Gamma_0r,Gamma_0a,...
    Gamma_1r,Gamma_1a,Gamma_ar,Gamma_aa,gamma_r,s0,s1,sr,sa,plus,minus,decay0r_1,decay1r_1,decayar_1,decay0a_1,decay1a_1,decayaa_1,...
    decay0r_2,decay1r_2,decayar_2,decay0a_2,decay1a_2,decayaa_2);

    for n=1:N
        pideal=reshape(rho(n,:),dim^2,dim^2);
        Population(n)= abs(trace(psi_0*psi_0'*pideal));
        Fidelity(n) = abs(trace(fidm*pideal));

        RHO1{jj,n} = pideal;
        population1{jj} = Population;
        fidelity1{jj} = Fidelity;
    end
    dm1 = pideal;
    % figure
    % legend
    hold on
    plot((jj-1)+tlist/time,Population,'b','LineWidth',2);
    plot((jj-1)+tlist/time,Fidelity,'r','LineWidth',2);

end
toc
grid minor