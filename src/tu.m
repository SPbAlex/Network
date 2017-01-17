clc
clear all
close all
f0 = 1200; % Гц
Vmod = 600; % Бод
Vinf = 1200; % бит/с
T = 1/600; % период
q = 2^(Vinf/Vmod); % количество сигналов
Q = sqrt(q); % промежуточная переменная
nfig = 1; % переменная для построения графиков
T0 = 1/f0; % период несущей
Ns = 20;   % число отсчетов на периоде несущей
E = 1; % энергия сигнала
A = (2*E/T)^(1/2); % максимальное абсолютное значение (амплитуда)
dt = min(1/f0)/Ns;
t = 0:dt:T;
phi1 = sqrt(2/T)*cos(2 * pi * f0 * t);
phi2 = sqrt(2/T)*sin(2 * pi * f0 * t);
% вывод на гарфик
figure(nfig);
nfig = nfig+1;
plot(t,phi1,'r.-', t,phi2,'b.-','lineWidth', 2);
% title('Вид базисных фукнций');
legend('cos(2Pft)','sin(2Pft)');
xlabel('t');
ylabel('|A|');
grid on;
s = zeros(q,length(t)); % инициализация массива СИГНАЛОВ размером q
figure (nfig);
nfig=nfig+1;
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        s(i+1,:) = A*(1 - 2*i1/(Q-1))*phi1 + A*(1 - 2*i2/(Q-1))*phi2;
        subplot(2,2,i+1);
        plot(t,s(i+1,:),'r.-','lineWidth', 2)
        title(['s_{', int2str(i), '}(t)']);
        xlabel('t,sec');
        ylabel('|s(t)|');
        grid on;
    end 
end 
df = (1/T)/20; % шаги по оси
f = [max(0,(f0-4/T)) : df : (f0+4/T)];
figure (nfig);
nfig=nfig+1;
plot(0,0);
S = zeros(q, length(f)); %массив спектров
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        S(i+1,:) = A*T/2*(A*(1-(2*i1)/(Q-1))*sinc((f-f0)*T)+sinc((f+f0)*T))+(A*T)/(2*j)*(A*(1-(2*i2)/(Q-1))*sinc((f-f0)*T)+sinc((f+f0)*T));
        subplot(2,2,i+1)
        plot(f, abs(S(i+1,:)), 'r.-','lineWidth', 2)
        title(['|S_{', int2str(i), '}(f)|']);
        xlabel('f,Hz');
        ylabel('|S(f)|');
        grid on;
        hold on;
    end
end

% ***********РАСЧЕТ И ПОСТРОЕНИЕ СПЕКТРОВ ВСЕХ СИГНАЛОВ (ОГИБАЮЩЕЙ)********
df = (1/T)/20; % шаги по оси
f = [max(0,(f0-4/T)) : df : (f0+4/T)];
% массив спектров функций
figure (nfig);
nfig=nfig+1;
plot(0,0);

S = zeros(q, length(f)); %массив спектров
for i1=0:Q-1
    for i2=0:Q-1
        i=i1*Q+i2;
        S(i+1,:) = A*T/2*(A*(1-(2*i1)/(Q-1))*sinc((f-f0)*T)+sinc((f+f0)*T))+(A*T)/(2*j)*(A*(1-(2*i2)/(Q-1))*sinc((f-f0)*T)+sinc((f+f0)*T));
        subplot(4,4,i+1)
        plot(f, abs(S(i+1,:)), 'r.-','lineWidth', 2)
        title(['|S_{', int2str(i), '}(f)|']);
        xlabel('f,Hz');
        ylabel('|S(f)|');
        grid on;
        hold on;
% выводим точки для определения полосы пропускания - [300 ; 3400]
        plot (300, 0, 'bs','lineWidth', 4);
        plot (3400, 0, 'bs','lineWidth', 4);
    end
end
%title('Виды спектров всех сигналов');

figure (nfig);
nfig=nfig+1;
N=25; % количество сигналов
ind=(floor(rand(1,N)*q)); % мультииндекс (выбор случайного сигнала)
S_spec=zeros(size(f)); % массив спектра последовательности
for l=0:N-1
    i=ind(l+1);
    S_spec=S_spec + S(i+1,:).*exp(j*2*pi*f*l*T);% получаем спектр последовательности
end
% вывод полученной последовательности на графике
plot(f,abs(S_spec), 'r.-','lineWidth', 2);
title('Спектор случайной последовательности');
grid on;
xlabel('f,Hz');
ylabel('|S_spec(f)|');
legend('S_spec(f)');

S_spektor = zeros(1,length(f));
p=0; % счетчик для построения спектра с узкой полосой
for i =[2,3,2,3,2,3,2,3,2,3];
    S_spektor = S_spektor + S(i+1,:).*exp(j*2*pi*f*p*T);
    p=p+1;
end
figure(nfig);
nfig = nfig+1;
plot(f,abs(S_spektor),'r.-','lineWidth', 2);
title('Вид спектра с узкой полосой');
grid on;
xlabel('f,Hz');
ylabel('|S_spektor(f)|');
legend('S_spektor(f)');

Spoint = zeros(q,2);% массив сигнального созвездия
for i=0:q-1
    Spoint(i+1,1) = sum(s(i+1,:).*phi1)*dt;
    Spoint(i+1,2) = sum(s(i+1,:).*phi2)*dt;
end;
% вывод на графике
figure(nfig);
nfig=nfig+1;
hold on;
for i=0:q-1
    plot(Spoint(i+1,1), Spoint(i+1,2), 'bo','lineWidth', 3)
    grid on
end;
% для красоты
for i=0:1:80
    plot(-40+i,0, 'r.','lineWidth', 1) % линия центру горизонтальная
    plot(0,40-i, 'r.','lineWidth', 1) % линия центру вертикальная
end
title('Сигнальное созвездие ');
xlabel('t,sec');
ylabel('|Spoint(t)|');
axis ('square');
dbmax = 5; %берем значения сигнал-шум =5
SNRdB = 4:1:10;
for j = 1:length(SNRdB)
    SNR = 10^(SNRdB(j)/10);
    sigma1 = sqrt(1/(2*SNR)*sum(sum(s.^2))/q);  %приложениe 2 конспекта лекций
end
sigma2 = 300;
% начаинаем строить облака рассеивания
figure(nfig);
nfig=nfig+1;
hold on; %накапливает точки на одном графике
r1 = zeros(q,length(t));
r2 = zeros(q,length(t));
isp = 2000; %задаем количество точек
in=0; %для цикла пробежки по значениям
% шум
for in=0:isp
    i = floor(q*rand);%случайное значение сигнала
    ro = s(i+1,:)+sigma1*randn(1,length(t));%смешиваем сигнал с шумом
    r1 = sum(sum(ro.*phi1))*dt;% первая точка массива
    r2 = sum(sum(ro.*phi2))*dt;% вторая точка массива
    plot(r1,r2,'b.','lineWidth', 2); %выводим точки на график
    grid on
end;
%центральные точки созвездия
for i=0:q-1
    plot(Spoint(i+1,1), Spoint(i+1,2), 'ro','lineWidth', 2)
end;
for i=0:1:160
    plot(-80+i,0, 'r.','lineWidth', 1) % линия центру горизонтальная
    plot(0,80-i, 'r.','lineWidth', 1) % линия центру вертикальная
end
title('Облака рассеивания ');
xlabel('t,sec');
ylabel('|ro(t)|');
axis ('square');%подбор масштаба

%проектирование оптимального приемника
%строим Ре путем эксперимента

dbmax = 5; %взять с шагом по графику
SNRdB = 0:1:10; % идем по графику от 0 до 10 с шагом 1
SNR = 10.^(SNRdB/10);  % переходим к децибелам (нужны только для вывода на графике)
% рис. 14.1 конспекта лекций

PeExper = zeros(1, length(SNRdB)); % объявляем одномерный массив для занесения значений практической вероятности ошибки
%      figure(nfig);
%    nfig = nfig + 1;
%    subplot(1,2,1)
    
    
SNRdBMax = 20; % максимальное значение точек =10
% NerrMax = 5;
d = 2*A/(sqrt(q)-1); % расстояние между сигналами
d=0.019441; % для нагляности сами рассчитали
ik = [[16 12 8 4];[15 11 7 3];[14 10 6 2];[13 9 5 1]]; % объявили двумерный массив с подобранными сигналами 
SNRdBExp = 0:1:SNRdBMax; % задали с интервалом от 0 до 20 и шагом 1
SNRExp = 10.^(SNRdBExp/10); % перевели в децибелы для рассчета сигмы...
sigma = sqrt( sum(sum(s.^2))./(2*SNRExp*q) ); % получение сигмы

PeExp = zeros(1,length(SNRdB)); % объявляем массив для практической верояности ошибки

for j = 1:length(SNRdB) % открываем главный цикл для получения значений верояности ошибки
    nRun = 0; % счетчик для цикла поиска ошибок
    nErr = 0; % счетчик подсчета ошибок
    nErrMax = 50; % считаем 
    while nErr < nErrMax % цикл поиска ошибок
        i = floor(q*rand)+1; % берем случайное значение от 1 до 16
    	r = s(i, :)+sigma(j)*randn(1, length(t)); % сигнал с шумом
        r1 = sum((r.*phi1)*dt); % точка для оси х
        r2 = sum((r.*phi2)*dt); % точка для оси у
 % определяем в какую решающую область попала точка
 % начинаем искать куда попала с координаты (-2;-2) - верхний левый угол
 % координатного созвездия
        r1k = max(-2, floor(r1/d)); % ищем в каком месте на оси х
        r1k = min(r1k, 1); % ищем на каком месте на оси у
        r1k = r1k+3; % получаем координату
        r2k = max(-2, floor(r2/d));% ищем в каком месте на оси х
        r2k = min(r2k, 1); % ищем на каком месте на оси у
        r2k = r2k+3; % получаем координату
 % проверяем на непопадание в решающую область
 % смотрим на значение случайно сгенерированного сигнала и смотрим на
 % координаты полученных значений 
        if (i ~= ik(r2k, r1k)) 
            nErr = nErr+1;
        end
 % счетчик проведенного теста
        nRun = nRun+1;
    end % закрываем поиск ошибочных точек
 % проверка полученных значений в выводе
    SNRdBExp(j); % значения сигнала шум
    nErr; % счетчика ошибок
    nRun; % счетчика проведенных тестов
    PeExp(j) = nErr/nRun; % получаем практическое значение вероятности ошибки
end % выход из главного цикла
    
% выводим график значений верояности ошибки
figure; % открываем новое окно
hold on; % делаем сеточку

%расчет теоретической вероятности ошибки
PeTheor = 1-(1-2*Qfun( sqrt( (3*SNR)/(q-1) ) )).^2;

% выводим в логирифмическом масштабе - используем функцию "semilogy"
semilogy(SNRdB, PeTheor, 'b', 'LineWidth', 2);
semilogy(SNRdB, PeExp, 'ro-', 'LineWidth', 2);
  

%semilogy(SNRdB, Pe, 'bo:')



