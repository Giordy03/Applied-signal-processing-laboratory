%% Exercise 1

clc;
clear;
close all;

f1 = 5;
f2 = 250;
B_D = 10;
fs = 1e3;
ts = 1/fs;
T_tot = 20;
N = T_tot/ts;
t = (0:N-1)*ts;

D3 = diric(2*pi*B_D*t, 3);
variance = 10;
dev = sqrt(variance);

W = dev*randn(1, N);
x = cos(2*pi*f1*t) + 3*(D3.^2).*cos(2*pi*f2*t) + W;

M = 25;
NFFT = N/M;
n_overlap = NFFT/2;
window = hamming(NFFT);
[Sx, f_x] = pwelch(x, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_x, 10*log10(abs(Sx)))
grid on
xlabel('Frequency (Hz)')
ylabel('10log_{10}|S_x (f)|')
title('PSD of Input Signal X(t)')

% create the filter
Wp = [210, 290]./(fs/2);
Ws = [190, 310]./(fs/2);
Rp = 1; %dB
Rs = 60; %attenuation in dB

% Butterworth filter

[n_order, Wn] = buttord(Wp,Ws,Rp,Rs);
cutoff = Wn*(fs/2);
[b, a] = butter(n_order, Wn, 'bandpass');

figure
freqz(b, a, 1024, fs)  
grid on
ylim([-100, 10])            
xlabel('Frequency (Hz)')    
ylabel('Magnitude (dB)')    
title('Butterworth Band-Pass Filter - Frequency Response')                   

y1 = filter(b, a, x);
[Sy1, f_y1] = pwelch(y1, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_y1, 10*log10(abs(Sy1)))
grid on
xlabel('Frequency (Hz)')
ylabel('10log_{10}|S_y (f)|')
title('Butterworth Band-Pass Filter - PSD of Output Signal Y(t)')

% Elliptic filter

[n_order_ell, Wn_ell] = ellipord(Wp,Ws,Rp,Rs);
cutoff = Wn_ell*(fs/2);

[b, a] = ellip(n_order_ell, Rp, Rs, Wn_ell, 'bandpass');

figure
freqz(b, a, 1024, fs)  
grid on
ylim([-100, 10])            
xlabel('Frequency (Hz)')    
ylabel('Magnitude (dB)')    
title('Elliptic Band-Pass Filter - Frequency Response')                   

y2 = filter(b, a, x);
[Sy2, f_y2] = pwelch(y2, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_y2, 10*log10(abs(Sy2)))
grid on
xlabel('Frequency (Hz)')
ylabel('10log_{10}|S_y (f)|')
title('Elliptic Band-Pass Filter - PSD of Output Signal Y(t)')

% FIR filter

f_fir = [Ws(1) Wp Ws(2)]*(fs/2);    
a = [0 1 0];

a_fir = 1;

dev = [10^(-Rs/20) 1-10^(-Rp/20) 10^(-Rs/20)]; 

[n,fo,ao,w] = firpmord(f_fir,a,dev,fs);

b = firpm(n,fo,ao,w, 'bandpass');

figure
freqz(b, a, 1024, fs)  
grid on
ylim([-100, 10])            
xlabel('Frequency (Hz)')    
ylabel('Magnitude (dB)')    
title('Equiripple FIR Band-Pass Filter - Frequency Response')  

y3 = filter(b, a_fir, x);
[Sy3, f_y3] = pwelch(y3, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_y3, 10*log10(abs(Sy3)))
grid on
xlabel('Frequency (Hz)')
ylabel('10log_{10}|S_y (f)|')
title('Equiripple FIR Band-Pass Filter - PSD of Output Signal Y(t)')

% Chebyshev

Wp = 6/(fs/2);
Ws = 14/(fs/2);
Rp = 1; %dB
Rs = 60; %attenuation in dB

[n_order_cheb, Wn_cheb] = cheb1ord(Wp,Ws,Rp,Rs);
cutoff = Wn_cheb*(fs/2);

[b_chev, a_chev] = cheby1(n_order_cheb, Rp, Wn_cheb, 'low');

figure
freqz(b_chev, a_chev, 1024, fs)  
grid on
ylim([-100, 10])    
xlim([0 50])
xlabel('Frequency (Hz)')    
ylabel('Magnitude (dB)')    
title('Chebyshev Type I Low-Pass Filter - Frequency Response') 

% Chebyshev my_filter()
M = 100;
delta = zeros(1,M);
delta(1) = 1; 
 
filtered_delta = my_filter(b_chev, a_chev, delta);

figure
impz(b_chev, a_chev, M);
grid on
hold on
stem(0:M-1, filtered_delta);
xlabel('Samples (n)')    
ylabel('Amplitude')   
legend('my\_filter()', 'impz()')
title('Chebyshev Type I Low-Pass Filter - Impulse Response Comparison') 

% my_filter() X var = 1
T_tot = 4;
N = T_tot/ts;
t = (0:N-1)*ts;

variance = 1;
dev = sqrt(variance);
D3 = diric(2*pi*B_D*t, 3);

W = dev*randn(1, N);
x = cos(2*pi*f1*t) + 3*(D3.^2).*cos(2*pi*f2*t) + W;

Y_4a = my_filter(b_chev, a_chev, x);
Y_4b = filter(b_chev, a_chev, x);

figure
plot(t, Y_4a, 'r-', 'LineWidth', 1.5) 
hold on
plot(t, Y_4b, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Time (s)')    
ylabel('Amplitude')   
legend('filter()', 'my\_filter()')
title('Filtering Comparison with Noise Variance 1 - Output Signal Y(t)')

% my_filter() X var = 0.1
T_tot = 4;
N = T_tot/ts;
t = (0:N-1)*ts;

variance = 0.1;
dev = sqrt(variance);
D3 = diric(2*pi*B_D*t, 3);

W = dev*randn(1, N);
x = cos(2*pi*f1*t) + 3*(D3.^2).*cos(2*pi*f2*t) + W;

Y_4a = my_filter(b_chev, a_chev, x);
Y_4b = filter(b_chev, a_chev, x);

figure
plot(t, Y_4a, 'r-', 'LineWidth', 1.5) 
hold on
plot(t, Y_4b, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Time (s)')    
ylabel('Amplitude')   
legend('filter()', 'my\_filter()')
title('Filtering Comparison with Noise Variance 0.1 - Output Signal Y(t)')

% CHECK ERROR BOX
b_error = [1, 1.3, 0.49, -0.013, -0.029];
a_error = [1, -0.4326, -1.6656, 0.1253, 0.2877];

z = roots(b_error);
p = roots(a_error);

figure
zplane(z,p)

error_disp = my_filter(b, a_error, x(1:N-1));