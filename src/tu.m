clc
clear all
close all
f0 = 1200;
Vmod = 600;
Vinf = 1200;
T = 1/600;
q = 2^(Vinf/Vmod);
Q = sqrt(q);
nfig = 1;
T0 = 1/f0;
Ns = 20;
E = 1;
A = (2*E/T)^(1/2);
dt = min(1/f0)/Ns;
t = 0:dt:T;
phi1 = sqrt(2/T)*cos(2 * pi * f0 * t);
phi2 = sqrt(2/T)*sin(2 * pi * f0 * t);

figure(nfig);
nfig = nfig+1;
plot(t,phi1,'r.-', t,phi2,'b.-','lineWidth', 2);

legend('cos(2Pft)','sin(2Pft)');
xlabel('t');
ylabel('|A|');
grid on;
s = zeros(q,length(t));
figure (nfig);
nfig=nfig+1;
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        s(i+1,:) = A*(1 - 2*i1/(Q-1))*phi1
			+ A*(1 - 2*i2/(Q-1))*phi2;
        subplot(2,2,i+1);
        plot(t,s(i+1,:),'r.-','lineWidth', 2)
        title(['s_{', int2str(i), '}(t)']);
        xlabel('t,sec');
        ylabel('|s(t)|');
        grid on;
    end 
end 
df = (1/T)/20;
f = [max(0,(f0-4/T)) : df : (f0+4/T)];
figure (nfig);
nfig=nfig+1;
plot(0,0);
S = zeros(q, length(f));
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        S(i+1,:) = A*T/2*(A*(1-(2*i1)/(Q-1))
		*sinc((f-f0)*T)
		+sinc((f+f0)*T))
		+(A*T)/(2*j)
		*(A*(1-(2*i2)/(Q-1))
		*sinc((f-f0)*T)+sinc((f+f0)*T));
        subplot(2,2,i+1)
        plot(f, abs(S(i+1,:)), 'r.-','lineWidth', 2)
        title(['|S_{', int2str(i), '}(f)|']);
        xlabel('f,Hz');
        ylabel('|S(f)|');
        grid on;
        hold on;
    end
end


df = (1/T)/20; 
f = [max(0,(f0-4/T)) : df : (f0+4/T)];

figure (nfig);
nfig=nfig+1;
plot(0,0);

S = zeros(q, length(f));
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        S(i+1,:) = A*T/2*(A*(1-(2*i1)/(Q-1))
		*sinc((f-f0)*T)
		+sinc((f+f0)*T))
		+(A*T)/(2*j)*(A*(1-(2*i2)/(Q-1))
		*sinc((f-f0)*T)+sinc((f+f0)*T));
        subplot(4,4,i+1)
        plot(f, abs(S(i+1,:)), 'r.-','lineWidth', 2)
        title(['|S_{', int2str(i), '}(f)|']);
        xlabel('f,Hz');
        ylabel('|S(f)|');
        grid on;
        hold on;

        plot (300, 0, 'bs','lineWidth', 4);
        plot (3400, 0, 'bs','lineWidth', 4);
    end
end


figure (nfig);
nfig=nfig+1;
N=25;
ind=(floor(rand(1,N)*q));
S_spec=zeros(size(f));
for l=0:N-1
    i=ind(l+1);
    S_spec=S_spec + S(i+1,:).*exp(j*2*pi*f*l*T);
end

plot(f,abs(S_spec), 'r.-','lineWidth', 2);

grid on;
xlabel('f,Hz');
ylabel('|S_spec(f)|');
legend('S_spec(f)');

S_spektor = zeros(1,length(f));
p=0; 
for i =[2,3,2,3,2,3,2,3,2,3];
    S_spektor = S_spektor 
		+ S(i+1,:).*exp(j*2*pi*f*p*T);
    p=p+1;
end
figure(nfig);
nfig = nfig+1;
plot(f,abs(S_spektor),'r.-','lineWidth', 2);

grid on;
xlabel('f,Hz');
ylabel('|S_spektor(f)|');
legend('S_spektor(f)');

Spoint = zeros(q,2);
for i=0:q-1
    Spoint(i+1,1) = sum(s(i+1,:).*phi1)*dt;
    Spoint(i+1,2) = sum(s(i+1,:).*phi2)*dt;
end;

figure(nfig);
nfig=nfig+1;
hold on;
for i=0:q-1
    plot(Spoint(i+1,1), Spoint(i+1,2), 'bo','lineWidth', 3)
    grid on
end;

for i=0:1:80
    plot(-40+i,0, 'r.','lineWidth', 1) 
    plot(0,40-i, 'r.','lineWidth', 1) 
end

xlabel('t,sec');
ylabel('|Spoint(t)|');
axis ('square');
dbmax = 5; 
SNRdB = 4:1:10;
for j = 1:length(SNRdB)
    SNR = 10^(SNRdB(j)/10);
    sigma1 = sqrt(1/(2*SNR)*sum(sum(s.^2))/q);  
end
sigma2 = 300;
figure(nfig);
nfig=nfig+1;
hold on; 
r1 = zeros(q,length(t));
r2 = zeros(q,length(t));
isp = 2000;
in=0; 

for in=0:isp
    i = floor(q*rand);
    ro = s(i+1,:)+sigma1*randn(1,length(t));
    r1 = sum(sum(ro.*phi1))*dt;
    r2 = sum(sum(ro.*phi2))*dt;
    plot(r1,r2,'b.','lineWidth', 2);
    grid on
end;

for i=0:q-1
    plot(Spoint(i+1,1), Spoint(i+1,2), 'ro','lineWidth', 2)
end;
for i=0:1:160
    plot(-80+i,0, 'r.','lineWidth', 1)
    plot(0,80-i, 'r.','lineWidth', 1)
end

xlabel('t,sec');
ylabel('|ro(t)|');
axis ('square');

dbmax = 5;
SNRdB = 0:1:10;
SNR = 10.^(SNRdB/10);

PeExper = zeros(1, length(SNRdB));
%      figure(nfig);
%    nfig = nfig + 1;
%    subplot(1,2,1)
    
    
SNRdBMax = 20; 
% NerrMax = 5;
d = 2*A/(sqrt(q)-1);
d=0.019441;
ik = [[16 12 8 4];[15 11 7 3];[14 10 6 2];[13 9 5 1]];
SNRdBExp = 0:1:SNRdBMax;
SNRExp = 10.^(SNRdBExp/10);
sigma = sqrt( sum(sum(s.^2))./(2*SNRExp*q) );

PeExp = zeros(1,length(SNRdB));

for j = 1:length(SNRdB)
    nRun = 0;
    nErr = 0;
    nErrMax = 50;
    while nErr < nErrMax
        i = floor(q*rand)+1;
    	r = s(i, :)+sigma(j)*randn(1, length(t)); 
        r1 = sum((r.*phi1)*dt); 
        r2 = sum((r.*phi2)*dt); 
        r1k = max(-2, floor(r1/d));
        r1k = min(r1k, 1);
        r1k = r1k+3;
        r2k = max(-2, floor(r2/d));
        r2k = min(r2k, 1);
        r2k = r2k+3;
        if (i ~= ik(r2k, r1k)) 
            nErr = nErr+1;
        end
        nRun = nRun+1;
    end

    SNRdBExp(j);
    nErr;
    nRun;
    PeExp(j) = nErr/nRun;
end
    
figure;
hold on;

PeTheor = 1-(1-2*Qfun( sqrt( (3*SNR)/(q-1) ) )).^2;

semilogy(SNRdB, PeTheor, 'b', 'LineWidth', 2);
semilogy(SNRdB, PeExp, 'ro-', 'LineWidth', 2);
