%% EX 3
clc, clear, close all

filename = 'pulse.txt';
delimiterIn = ' ';
headerlinesIn = 1;
Data_struct = importdata(filename,delimiterIn,headerlinesIn);
Led_R = Data_struct.data(:,1); % RED (R)
Led_IR = Data_struct.data(:,2); % INFRARED (IR)

fs = 100;
ts = 1/fs;

t_start = 10;
N_start = t_start*fs;
t_end = 60;
N_end = (t_end+t_start)*fs;
Led_R = Led_R(N_start:N_end);
Led_IR = Led_IR(N_start:N_end);

t = 0:ts:t_end;

figure
sgtitle('Read signal')
subplot(2,1,1)
plot(t, Led_R)
xlabel("time [s]"), ylabel("Amplitude")
title('RED')
grid on
subplot(2,1,2)
plot(t, Led_IR)
xlabel("time [s]"), ylabel("Amplitude")
title('IRED')
grid on

%% design IIR filter
omega_p = 3;
wp = omega_p*2*pi;
omega_s = 6;
ws = omega_s*2*pi;
Rp = 1; %dB
As = 60; %dB

N = ceil(log10((10^(Rp/10) - 1)/(10^(As/10) - 1))/(2*log10(wp/ws)));
fprintf("The Butterworth low-pass filter order is: %d\n", N)
wc_1 = wp/((10^(Rp/10)-1)^(1/(2*N)));
wc_2 = ws/((10^(As/10)-1)^(1/(2*N)));
wc = mean([wc_1, wc_2]);
omega_c = round(wc/(2*pi), 3);
fprintf("The Butterworth low-pass filter cut-off frequency is: %f Hz\n", omega_c)
f = linspace(0, 40, 1001);

H_squared = 1./(1 + (f./omega_c).^(2*N));

figure
plot(f, 10*log10(H_squared))
title("Squared magnitude of the analog filter")
xlabel("Frequency [Hz]"), ylabel("10log_{10}(|H_a(j\Omega)|^2) [dB]")
grid on

n = 0:2*N-1;
p_k = wc*exp(1j*n*pi/N)*exp(1j*(pi*(N+1)/(2*N)));

figure

plot(p_k, 'o')
grid on
axis('equal')
xlabel("Real part"), ylabel("Imaginary part")
title("All poles of the designed IIR filter")

p_k = p_k(real(p_k)<0);
%%
[num_Ha, den_Ha] = zp2tf([], p_k, wc^N);

figure
zplane(num_Ha, den_Ha)
grid on

[num_Hz, den_Hz] = impinvar(num_Ha, den_Ha, fs);

figure
freqz(num_Hz, den_Hz, 1024, fs);

R_filt = filtfilt(num_Hz, den_Hz, Led_R);
IR_filt = filtfilt(num_Hz, den_Hz, Led_IR);

%% create FIR filter using hamming window
omega_s = 0.05/fs;
omega_p = 0.75/fs;
Bt = abs(omega_p - omega_s)*2*pi;

N = ceil(6.6*pi/Bt);
fprintf('The number of points is %d.\n', N)

fc_norm = (omega_p + omega_s)/(2);
fprintf("The cut-off frequency of the high-pass FIR filter with Hamming window is %f Hz \n", fc_norm*fs)
n = 0:(N-1);
w_Hamming = 0.54 - 0.46*cos(2*pi*n/(N-1));

M = ceil((N-1)/2);

h_id = (n==M) - 2*fc_norm*sinc(2*fc_norm*(n-M));
h = h_id.*w_Hamming;

figure
freqz(h, 1, 1024, fs)

[H, f_freqz] = freqz(h, 1, 1024, fs);
figure
plot(f_freqz, 20*log10(abs(H)))
grid on
xlim([0, 3])
title("FIR high-pass enlargement")
xlabel("Frequency [Hz]"), ylabel("Magnitude [dB]")

R_filt2 = filtfilt(h, 1, R_filt); % perchè è un fir quindi non ha poli ed è tutto numeratore 
IR_filt2 = filtfilt(h, 1, IR_filt);

figure
sgtitle("RED signal")
subplot(3,1,1)
plot(t, Led_R)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("original")
subplot(3,1,2)
plot(t, R_filt)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("After low pass")
subplot(3,1,3)
plot(t, R_filt2)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("After low pass and high pass")

figure
sgtitle("INFRARED signal")
subplot(3,1,1)
plot(t, Led_IR)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("Original")
subplot(3,1,2)
plot(t, IR_filt)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("After low pass")
subplot(3,1,3)
plot(t, IR_filt2)
grid on
xlabel("time [s]"), ylabel("Amplitude")
title("After low pass and high pass")
% Compute FFT

N_fft = 2^nextpow2(length(R_filt2));
RED_fft = fftshift(fft(R_filt2, N_fft)); 

df = fs/N_fft;
f = -fs/2:df:fs/2-df;

[max_val, idx_max] = max(RED_fft);
f_max = f(idx_max);
BPM = f_max*60;

figure
plot(f, 20*log10(abs(RED_fft)))
grid on
title(sprintf("Filtered RED signal - BPM = %.3f", BPM))

%% step 7 - saturation computation
% for DC signal use the low pass filter output
[DC_red_val, DC_red_loc] = findpeaks(-R_filt);
I_red_DC = interp1(t(DC_red_loc), -DC_red_val, t, "spline");

[DC_infra_val, DC_infra_loc] = findpeaks(-IR_filt);
I_infra_DC = interp1(t(DC_infra_loc), -DC_infra_val, t, "spline");

% for AC signal use the high pass filter output
[AC_min_red_val, AC_min_red_loc] = findpeaks(-R_filt2);
[AC_max_red_val, AC_max_red_loc] = findpeaks(R_filt2);
I_red_AC = interp1([t(AC_min_red_loc), t(AC_max_red_loc)], [-AC_min_red_val; AC_max_red_val], t, "spline");
I_red = interp1([t(AC_min_red_loc), t(AC_max_red_loc)], [-DC_red_val; -DC_red_val+ AC_min_red_val + AC_max_red_val], t, "spline");

[AC_min_infra_val, AC_min_infra_loc] = findpeaks(-IR_filt2);
[AC_max_infra_val, AC_max_infra_loc] = findpeaks(IR_filt2);
I_infra_AC = interp1([t(AC_min_infra_loc), t(AC_max_infra_loc)], [-AC_min_infra_val; AC_max_infra_val], t, "spline");
I_infra = interp1([t(AC_min_infra_loc), t(AC_max_infra_loc)], [-DC_infra_val; -DC_infra_val+ AC_min_infra_val + AC_max_infra_val], t, "spline");

figure
sgtitle("Peaks")
subplot(2, 2, 1)
hold on
plot(t(DC_red_loc), -DC_red_val, "xk")
title("DC red minimum")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 2)
hold on
plot(t(AC_min_red_loc), -AC_min_red_val, "xk")
plot(t(AC_max_red_loc), AC_max_red_val, "or")
title("AC red extrema")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 3)
hold on
plot(t(DC_infra_loc), -DC_infra_val, "or")
title("DC infrared extrema")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 4)
hold on
plot(t(AC_min_infra_loc), -AC_min_infra_val, "xk")
plot(t(AC_max_infra_loc), AC_max_infra_val, "or")
title("AC infrared extrema")
xlabel("Time [s]"); ylabel("Amplitude")



deltaN = 1*fs;
n = 1:deltaN:length(I_red_DC);

R = (I_red_AC(n)./I_red_DC(n)./(I_infra_AC(n)./I_infra_DC(n)));
sat = 110 - 25*mean(R);

figure
sgtitle(sprintf("SaO_2 = %.4f%%", sat))
subplot(2, 2, 1)
hold on
plot(t(DC_red_loc), -DC_red_val, "or")
plot(t, I_red, "-k")
ylim([2e5, 2.02e5])
title("DC red signal")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 2)
hold on
plot(t, I_red_AC, "-k")
plot(t(AC_min_red_loc), -AC_min_red_val, "or")
plot(t(AC_max_red_loc), AC_max_red_val, "or")
ylim([-400, 300])
title("AC red signal")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 3)
hold on
plot(t(DC_infra_loc), -DC_infra_val, "or")
plot(t, I_infra, "-k")
ylim([2.62e5, 2.65e5])
title("DC infrared signal")
xlabel("Time [s]"); ylabel("Amplitude")
subplot(2, 2, 4)
hold on
plot(t, I_infra_AC, "-k")
plot(t(AC_min_infra_loc), -AC_min_infra_val, "or")
plot(t(AC_max_infra_loc), AC_max_infra_val, "or")
ylim([-600, 700])
title("AC infrared signal")
xlabel("Time [s]"); ylabel("Amplitude")
