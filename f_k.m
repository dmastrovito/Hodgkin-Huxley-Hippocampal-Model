clear all;
cput = cputime;

N = 10; % typically 100 number of cells in network
Msyn= 10; % experimentally determined fixed average number of synaptic inputs per neuron
p = Msyn/N; % probability of two neurons being connected
Esyn = -75;
ps=random('bino',1,.6,N,N);
ps = ones(N,N);

alpha = 12; % msec^-1
beta = .1; %msec^-1
s_inf = alpha/(alpha + beta);

%Iapp from a gaussian distribution with mean u and std sigma -- so wait
%this is periodic input with freq 35Hz
Iu = 1;
Isigma = 0:.03/10:.03;
%Iapp = .5:2:20;

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

f = zeros(N,1);
fsigma=zeros(length(Isigma),1);
fu = zeros(length(Isigma),1);
K = zeros(length(Isigma),1);

for Is=1:length(gsyn)
    spiketrain = zeros(N,T/dt);
    
    for t=2:T/dt
    
        V_last = V(:,t-1); 
        s_last = s(:,t-1);
    
        for i=1:N
        %if (t > 1000/.5)
            I = Iapp;
        %end
            [V(i,t),n(i,t),h(i,t),s(i,t)] = hodgkin_huxley_synapse(dt,[n(i,t-1),h(i,t-1)],V_last,s_last,i,I,Esyn,ps(i,:),gsyn(Is));
            if ((V(i,t) > 0) && (V(i,t-1) < 0)) 
                spiketrain(i,t) = 1;
            end
        end
        
    end
    
    
    Knum=zeros(0,N);
    Tau = 2/dt; % 1msec/dt

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

K(Is) = mean2(not(isnan(Knum)))
    
    
end



colors = ['k'; 'b' ;'g'; 'y'; 'r'];

figure

plot(cm/gsyn, K)
xlabel('Tau syn')
ylabel('Coherence')

