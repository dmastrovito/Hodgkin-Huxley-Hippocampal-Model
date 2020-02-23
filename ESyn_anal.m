clear all;
cput = cputime;

N = 50; % typically 100 number of cells in network
Msyn= 60; % experimentally determined fixed average number of synaptic inputs per neuron
p = Msyn/N; % probability of two neurons being connected
ps = ones(N,N);
Esyn = [-90, -80, -60,-50,-30,-10];

alpha = 12; % msec^-1
beta = .1; %msec^-1
s_inf = alpha/(alpha + beta);

%Iapp from a gaussian distribution with mean u and std sigma 
Iu = 1;
Isigma = .03;
Iapp = normrnd(Iu,Isigma,N,1);
Iapp=1.2;

dt = .05; %ms
T = 2000;

V = zeros(N,T/dt);
V_last = zeros(N);

n = zeros(N,T/dt);
h = zeros(N,T/dt);
s = zeros(N,T/dt);
s_last = zeros(N);



%Setting Initial Conditions
%uniform distribution of V(1) between -50 and -70
%V(:,1) = -70 + (-50-(-70)).*rand(10,1);
V(:,1) = -1*randi([50,70],N,1);

alpha_h  = .07.*exp(-(V(:,1)+58)/20); 
beta_h = 1/(exp(-.1.*(V(:,1)+28)) +1);
h_inf  = alpha_h'./(alpha_h' + beta_h);
h(:,1) = h_inf;


alpha_n = -.01.*(V(:,1) +34)./(exp(-.1.*(V(:,1)+34)) -1);
beta_n = .125.*exp(-(V(:,1)+44)/80);
n_inf = alpha_n./(alpha_n + beta_n);
n(:,1) =  n_inf;
s(:,1) = s_inf;


%Build network connection matrix
%randint(N,1,[0 1])
K = zeros(length(Esyn),1);
I=Iapp;

for e=1:length(Esyn)
    spiketrain = zeros(N,T/dt);
    for t=2:T/dt
    
        V_last = V(:,t-1); 
        s_last = s(:,t-1);
    
        for i=1:N
        %if (t > 1000/.5)
            I = Iapp;
        %end
            [V(i,t),n(i,t),h(i,t),s(i,t)] = hodgkin_huxley_synapse(dt,[n(i,t-1),h(i,t-1)],V_last,s_last,i,I,Esyn(e),ps(i,:));
            if ((V(i,t) > 0) && (V(i,t-1) < 0)) 
                spiketrain(i,t) = 1;
            end
        end
    end



time=1:T/dt;
colors = ['k'; 'b' ;'g'; 'y'; 'r'];
if N >5
    colors =repmat(colors, N/5,1);
end

figure
for i=1:N
    %subplot(5,2,i)
    plot(time(1:10000)*dt,V(i,1:10000),colors(i))
    hold on
end
xlabel('time msec')
ylabel('Voltage mV')
title(strcat('Network Activity ',num2str(N), ' Neurons'))
legend(strcat("Esyn = ",num2str(Esyn(e))))



Knum=zeros(0,N);
Tau = 1/dt; % 1msec/dt

for i=1:N
    for j=1:N
            ith = zeros(T/dt/2,1);
            jth = zeros(T/dt/2,1);
            for tau=20000:Tau:T/dt-1 %tau = 1 msec
                ith(tau) = sum(spiketrain(i,tau:tau+1));
                jth(tau) = sum(spiketrain(j,tau:tau+1));
            end
            Knum(i,j)= (sum(ith.*jth))/sqrt(sum(ith)*sum(jth));
    end
end

K(e) = mean2(Knum)
end
cput = cputime - cput
figure
plot(Esyn,K);
xlabel('Es in mV')
ylabel('Network Coherence K')
title(strcat('Network Synchrony as a function of Esyn  Run Time=',num2str(cput),' seconds'))

