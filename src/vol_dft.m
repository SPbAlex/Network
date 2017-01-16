    N = 8;
    n = 0:N-1;
    k = 1;

    %%%% FIRST TASK: 1 freq %%%%
    disp('FIRST TASK: 1 freq');
    signal = cos(2*pi*n*k/N) + 1i*sin(2*pi*n*k/N);
    plot(n, signal);
    disp('FFT result');
    fft_res = fft(signal)
    %%%% END TASK %%%%

    %%%% SECOND TASK: 2 freq %%%%
    disp('SECOND TASK: 2 freq');
    signal1 = cos(2*pi*n*2/N) + 1i*sin(2*pi*n*2/N);
    signal2 = cos(2*pi*n*5/N) + 1i*sin(2*pi*n*5/N);
    signal = signal1 + signal2;
    plot(n, signal);
    disp('FFT result');
    fft_res = fft(signal)
    %%%% END TASK %%%%

    %%%% THIRD TASK: ifft -> sin() %%%%
    disp('THIRD TASK: ifft sin()');
    signal = sin(2*pi*n/N);
    signal = [0 + 8* 1i, 0 - 0*1i,0 - 0*1i,0 - 0*1i,0 - 0*1i, ...
              0 - 0*1i,0 - 0*1i,0 + 8*1i]
    plot(n, signal);
    disp('IFFT result');
    ifft_res = ifft(signal)
    plot(n, ifft_res);
    %%%% END TASK %%%%