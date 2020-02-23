function [V,n,h,s] = Hodgkin_huxley_synapse(dt,probs,Vs,Ss,i,Iapp,Esyn,ps,Msyn)

Cm = 1; %uF/cm^2

gNa = 35; %mS/cm^2
gL = .1; %mS/cm^2
gK = 9; %mS/cm^2
gsyn = .1; %mS/cm^2

ENa = 55; %mV
EL  = -65; % mV;
EK = -90; %mV
%Esyn = -75; %mV inhibitory
%Esyn = 0; % mV excitatory

n = probs(1);
h = probs(2);
s = Ss(i);
V = Vs(i);

phi = 5;
alpha = 12; % msec^-1
beta = .1; %msec^-1

alpha_m = -.1*(V+ 35)/(exp(-.1*(V+35))-1);
beta_m = 4*exp(-(V+60)/18);
m_inf = alpha_m/(alpha_m + beta_m);


alpha_h  = .07*exp(-(V+58)/20);
beta_h = 1/(exp(-.1*(V+28)) +1);
h1 = phi*(alpha_h*(1-h) - beta_h*h);
h1h = h + h1*dt;
h2 = phi*(alpha_h*(1-h1h) - beta_h*h1h);
h = h + (h1 + h2)/2*dt;

alpha_n = -.01*(V +34)/(exp(-.1*(V+34)) -1);
beta_n = .125*exp(-(V+44)/80);
n1 = phi*(alpha_n*(1-n)-beta_n*n);
n1n = n + n1*dt;
n2 = phi*(alpha_n*(1-n1n)-beta_n*n1n);
n = n + (n1 +n2)/2*dt;

INa = gNa*m_inf^3*h*(V-ENa);
IL = gL*(V - EL);
IK = gK*n^4*(V-EK);


theta = 0; %mV
F = 1/(1 + exp(-(V - theta)/2));
s1 = (alpha * F*(1-s)) - beta*s;
s1s = s + s1*dt;
s2 = (alpha * F*(1-s1s)) - beta*s1s;
s = s + (s1 + s2)/2*dt;

Ss(i) = s;

Isyn = sum(ps'.*gsyn/Msyn.*Ss.*(Vs-Esyn));

V1 = (-INa - IK - IL - Isyn + Iapp)/Cm;
V1V = V + V1*dt;

INa = gNa*m_inf^3*h*(V1V-ENa);
IL = gL*(V1V - EL);
IK = gK*n^4*(V1V-EK);
Vs(i) = V1V;
Isyn = sum(ps'.*gsyn/Msyn.*Ss.*(Vs-Esyn));
V2 = (-INa - IK - IL - Isyn + Iapp)/Cm;

V = V + (V1 + V2)/2*dt;



end
