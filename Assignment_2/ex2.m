%% EXERCISE 2

clc; clear all;
close all;

%% Bandpass

fprintf('Bandpass: Hamming\n')
fs = 1e3;
%from table -> min attenuation: 50 dB -> Hamming
trans_band = 30;
trans_band_norm = trans_band/(fs)*2*pi;
N = ceil(6.6*pi/trans_band_norm);
fprintf('The number of points is %d.\n', N)

f_cent = 200;
passband_length = 70;

Wp = [f_cent - passband_length/2, f_cent + passband_length/2]/(fs)*2*pi;
Ws = [f_cent - passband_length/2 - trans_band, f_cent + passband_length/2 + trans_band]/(fs)*2*pi;
fc_norm = (Wp + Ws)/(2*2*pi);

fprintf('The normalized (with respect to f_s) cutoff frequencies are %.2f and %.2f.\n\n', fc_norm)

n = 0:(N-1);
w_Hamming = 0.54 - 0.46*cos(2*pi*n/(N-1));
df = fs/N;
f = -fs/2:df:fs/2-df;

M = ceil((N-1)/2);

h_id = 2*fc_norm(2)*sinc(2*fc_norm(2)*(n-M)) - 2*fc_norm(1)*sinc(2*fc_norm(1)*(n-M)); 
h = h_id.*w_Hamming;

imp_res = impz(h,1,1024, fs);
[H, f] = freqz(imp_res, 1, 1024, fs);
n_imp = 0:length(imp_res)-1;

figure
plot(n_imp, imp_res)
grid on
xlabel('n (samples)')
ylabel('h[n]')
title('Impulse response of the Hamming bandpass filter')
%saveas(gcf, 'dataes2/impres.png');

figure
plot(f, 20*log10(abs(H)))
grid on
hold on
xline(fc_norm(1)*fs, '--r', 'LineWidth', 1.2)
xline(fc_norm(2)*fs, '--r', 'LineWidth', 1.2)
xline(165, ':k', 'LineWidth', 1.2)
xline(235, ':k', 'LineWidth', 1.2)
xline(135, ':k', 'LineWidth', 1.2)
xline(265, ':k', 'LineWidth', 1.2)
xlabel('f [Hz]')
ylabel('|H(f)|^2_{dB}')
title('Frequency response of the Hamming bandpass filter')
%saveas(gcf, 'dataes2/freqres.png');



%% Bandstop

fprintf('Bandstop: Blackman\n')
%from table -> min attenuation: 70 dB -> Blackman
trans_band = 40;
trans_band_norm = trans_band/(fs)*2*pi;
N = ceil(11*pi/trans_band_norm);
fprintf('The number of points is %d.\n', N)

Wp = [0, 255, 395, 500]/(fs)*2*pi;
Ws = [295, 355]/(fs)*2*pi;
fc_norm = (Wp(2:3) + Ws)/(2*2*pi);

fprintf('The normalized (with respect to f_s) cutoff frequencies are %.2f and %.2f.\n', fc_norm)

n = 0:(N-1);
w_Blackman = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1));

df = fs/N;
f = -fs/2:df:fs/2-df;
f_cent_norm = 250/(fs/2);

M = ceil((N-1)/2);

h_id = (n==M) - (2*fc_norm(2)*sinc(2*fc_norm(2)*(n-M)) - 2*fc_norm(1)*sinc(2*fc_norm(1)*(n-M))); 
h = h_id.*w_Blackman;

imp_res = impz(h,1, 1024, fs);
[H, f] = freqz(imp_res, 1, 1024, fs);
n_imp = 0:length(imp_res)-1;

figure
plot(n_imp, imp_res)
grid on
xlabel('n (samples)')
ylabel('h[n]')
title('Impulse response of the Blackman bandstop filter')
%saveas(gcf, 'dataes2/impresstop.png');


figure
plot(f, 20*log10(abs(H)))
grid on
hold on
xline(fc_norm(1)*fs, '--r', 'LineWidth', 1.2)
xline(fc_norm(2)*fs, '--r', 'LineWidth', 1.2)
xline(255, ':k', 'LineWidth', 1.2)
xline(395, ':k', 'LineWidth', 1.2)
xline(355, ':k', 'LineWidth', 1.2)
xline(295, ':k', 'LineWidth', 1.2)
xlabel('f [Hz]')
ylabel('|H(f)|^2_{dB}')
title('Frequency response of the Blackman bandstop filter')
%saveas(gcf, 'dataes2/freqresstop.png');

%% Kaiser

% HIGH-PASS KAISER
As = 40;
Bt = 200;
fc = 1.6e3;
fs = 4e3;
Wp = (fc+Bt/2)/(fs)*2*pi;
Ws = (fc-Bt/2)/(fs)*2*pi;

beta = 0.5842*(As-21)^0.4 + 0.07886*(As-21);

N = ceil((As-8)/(2.285*(abs(Wp-Ws))));

[Bk_1, N, beta]  = my_Kaiser_filter(As, Bt, fs, fc, '-hp');
b1= fir1(N-1, fc/(fs/2), kaiser(N, beta), 'high');

[H1, ~] = freqz(Bk_1, 1, 1024, fs);

[H2, f] = freqz(b1, 1, 1024, fs);


figure
plot(f, 20*log10(abs(H1)))
grid on
hold on
plot(f, 20*log10(abs(H2)), "--r")
xlabel('f [Hz]')
ylabel('|H(f)|^2_{dB}')
legend('my\_Kaiser\_filter', 'fir1', 'Location', 'northwest')
title('Comparison between frequency responses')
%saveas(gcf, 'dataes2/comphigh.png');


% BAND-PASS KAISER
As = 60;
Bt = 100;
fc = [800, 1.2e3];

Wp = (fc+[+Bt/2, -Bt/2])/(fs)*2*pi;
Ws = ([0, fc(1)-Bt/2, fc(2)+Bt/2, fs/2])/(fs)*2*pi;

beta = 0.1102*(As-8.7);

N = ceil((As-8)/(2.285*(abs(Wp(1)-Ws(2)))));

[Bk_1, N_1, beta_1]  = my_Kaiser_filter(As, Bt, fs, fc, '-bp');
b1=fir1(N-1, fc/(fs/2), kaiser(N, beta), 'bandpass');

[H1, ~] = freqz(Bk_1, 1, 1024, fs);

[H2, f] = freqz(b1, 1, 1024, fs);


figure
plot(f, 20*log10(abs(H1)))
grid on
hold on
plot(f, 20*log10(abs(H2)), "--r")
xlabel('f [Hz]')
ylabel('|H(f)|^2_{dB}')
legend('my\_Kaiser\_filter', 'fir1')
title('Comparison between frequency responses')
%saveas(gcf, 'dataes2/compband.png');
