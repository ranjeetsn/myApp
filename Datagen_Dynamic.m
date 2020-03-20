clc
clear all
close all
rng(1)

%% Estimated model and error-variances obtained using DIPCA in state-space form
%[A,B,C,D]=tf2ss([-0.0186 1.1490 1.6676 0.0927 -0.0639 -0.0378],[1 0.1965 0.0294 -0.0304 0.0188 0.5744]); %correct order=5
%D=0;
%sigma=[0.4351 0.1214]; %correct order=5


%% Simulating the same data used to obtain the above model from DIPCA - can also do this by loading .mat file into the workspace
N=1023;
u1=idinput(N,'prbs');
y1=zeros(N,1);

% Simulating the data using DGP
order=2;
for k=(order+1):N
       y1(k)=-0.2*y1(k-1)-0.6*y1(k-2)+1.2*u1(k-1)+1.6*u1(k-2);
end

ek1=0.280*wgn(N,1,1); %error in u
ek2=0.638*wgn(N,1,1); %error in y

u=u1+ek1; y=y1+ek2;
snr=[var(u1)/var(ek1) var(y1)/var(ek2)];
%z1 = [u1 y1];
z1 = [y1 u1];
z  = [y u];