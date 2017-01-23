%-------------------------------------------------------------------------%
% clear all
% close all
clc

q = 2;
K = 200;
SNRdB = 0:6;
SNR = 10.^(SNRdB/10);

PeExp = zeros(1,length(SNRdB));
SNR_exp_dB = zeros(1,length(SNRdB));
SNR_exp = zeros(1,length(SNRdB));
N = zeros(1,length(SNRdB));
s = zeros(q,K);

s(1,:) = 1;
s(2,:) = -1;
for i = 1:length(SNRdB)
    step = 50000;
    nError = 0;
    SNR_sum = 0;
    for k = 1:step
        
        M = randi(2)-1;

        if M == 1
            s_message = s(1,:);
        else
            s_message = s(2,:);
        end

        tmp = sum(sum(s.^2));
        sigma = sqrt(tmp / (2 * SNR(i) * q));
        noise = sigma * randn(1,K);

        R = s_message + noise;

        if sum(R) >= 0
            Y = 1;
        else
            Y = 0;
        end

        if Y ~= M
            nError = nError + 1;
        end
        E = sum(s_message.^2);
        Enoise = sum(noise.^2);
        SNR_step = E / Enoise;
        SNR_sum = SNR_sum + SNR_step; 
    end
    disp(SNRdB(i));
    SNR_exp(i) = (K / q) * SNR_sum / step;
%     SNR_exp_dB(i) = 10 * log10(SNR_exp(i));
    PeExp(i) = nError / step;
    end
figure(1);
hold on
plot(SNRdB,PeExp,'r','linewidth',2);
grid on
legend('PeExp')
figure(2);
hold on
plot(SNR,SNR_exp,'b','linewidth',2);
grid on
legend('SNR exp')