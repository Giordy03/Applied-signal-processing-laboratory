clc
clear
close all
format short e

% gaussian filter
%n_overlap -> absolute parameter: not in percentage
% we need to be able to distinguish the two close deltas of the cosine

% simple periodogram -> no use other command, only pwelsh with set
% parameters St = n, solo un segmento

%BARTLETT: nfft adjusted to be equal to D

% my_DFT -> create a circulant matrix and then perform the matrix by vector
% operation -> no for loops

%keep the functions local, NOT in another file

nu = 0.01;
f_1 = 100;
f_2 = 500;
f_3 = 510;

fc = 2e3;
T_tot = 20;

tc = 1/fc;
N = T_tot*fc;
t = (0:N-1)*tc;

dev = sqrt(25);
w = dev*randn(1,N);

x = 20/(sqrt(2*pi)*nu).*exp(-(t.^2)./(2*nu^2)).*cos(2*pi*f_1*t)+cos(2*pi*f_2*t)+cos(2*pi*f_3*t)+ w;

NFFT = 128;
n_overlap = 0;
window = hamming(NFFT);
[Sx, f] = pwelch(x, window, n_overlap, NFFT, fc, 'centered');

figure
plot(f, 10*log10(abs(Sx)))
ylabel('10log_{10}|S_x(f)|')
xlabel('Frequency (Hz)')
title('PSD in logarithmic scale')

var_Sx = var(Sx);

fprintf('The estimation variance, using a Hamming window, NFFT = %d, and n_overlap = %d is: %d \n', NFFT, n_overlap, var_Sx)

NFFT = 1024;
n_overlap = NFFT/2;
window = hamming(NFFT);
[Sx, f] = pwelch(x, window, n_overlap, NFFT, fc, 'centered');

figure
plot(f, 10*log10(abs(Sx)))
ylabel('10log_{10}|S_x(f)|')
xlabel('Frequency (Hz)')
title('PSD in logarithmic scale')

var_Sx = var(Sx);

fprintf('The estimation variance, using a Hamming window, NFFT = %d, and n_overlap = %d is: %d \n', NFFT, n_overlap, var_Sx)

% simple periodogram
NFFT = N;
n_overlap = 0;
window = rectwin(NFFT);
[Sx_simple, f_simple] = pwelch(x, window, n_overlap, NFFT, fc, 'centered');

% Bartlett
M = 25;
NFFT = floor(N/M);
n_overlap = 0;
window = rectwin(NFFT);
[Sx_bartlett, f_bartlett] = pwelch(x, window, n_overlap, NFFT, fc, 'centered');

% Welch
NFFT = floor(N/M);
n_overlap = NFFT/2;
window = hann(NFFT);
[Sx_welch, f_welch] = pwelch(x, window, n_overlap, NFFT, fc, 'centered');

figure
plot(f_simple, 10*log10(abs(Sx_simple)))
hold on
grid on
plot(f_bartlett, 10*log10(abs(Sx_bartlett)))
plot(f_welch, 10*log10(abs(Sx_welch)))
legend('Simple', 'Bartlett', 'Welch')
ylabel('10log_{10}|S_x(f)|')
xlabel('Frequency (Hz)')
title('Power Spectral Density: Simple vs Bartlett vs Welch Periodogram');


fprintf(['The variance of Simple periodogram is ', num2str(var(Sx_simple)), '\n'])
fprintf(['The variance of Bartlett periodogram is ', num2str(var(Sx_bartlett)), '\n'])
fprintf(['The variance of Welch periodogram is ', num2str(var(Sx_welch)), '\n'])


[R_N, lags] = xcorr(x, 'unbiased');
[R_prime_N, lags_prime] = xcorr(x, 'biased');

figure
subplot(2,1,1)
plot(lags, 10*log10(abs(R_N)));
xlabel('Lag')
ylabel('10log_{10}|R_N|')
title('Unbiased Autocorrelation Estimator')

subplot(2,1,2)
plot(lags_prime, 10*log10(abs(R_prime_N)));
xlabel('Lag')
ylabel('10log_{10}|R''_N|')
title('Biased Autocorrelation Estimator')

N = 165;
M = floor(N/2);
n = (-M:M)';
h = 0.3*sinc(0.3*n).*hamming(N);
y = filter(h, 1, x);

% welch on y
D = 1600;
n_overlap = D/2;
window = rectwin(D);
Sy_rect = pwelch(y, window, n_overlap, D, fc, 'centered');
window = hamming(D);
Sy_hamm = pwelch(y, window, n_overlap, D, fc, 'centered');
window = hann(D);
[Sy_hann, f] = pwelch(y, window, n_overlap, D, fc, 'centered');

figure;
plot(f, 10*log10(abs(Sy_rect)), 'LineWidth', 1.2); 
hold on
grid on
plot(f, 10*log10(abs(Sy_hamm)), 'LineWidth', 1.2);
plot(f, 10*log10(abs(Sy_hann)), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('10log_{10}|S_y(f)|');
title('Welch Periodogram of Filtered Signal with Different Windows');
legend('Rectangular Window', 'Hamming Window', 'Hann Window', 'Location', 'Best');

M_bartlett = 25;
[Sy_myBartlett, f1] = my_PSD(y, length(y), M_bartlett);
[Sy_mySimple, f2] = my_PSD(y, length(y));

SyM_myBartlett = mean(Sy_myBartlett, 1);

window = rectwin(length(y));
Sy_simple = pwelch(y, window, 0, length(y), 1, "centered");
window = rectwin(length(y) / M_bartlett);
Sy_Bartlett = pwelch(y, window, 0, length(y) / M_bartlett, 1, "centered");

figure
plot(f1, 10*log10(abs(Sy_Bartlett)), 'LineWidth', 1.2)
hold on
grid on
plot(f1, 10*log10(abs(SyM_myBartlett)), 'LineWidth', 1.2)
xlabel('Frequency (Hz)');
ylabel('10log_{10}|S_y(f)|');
title('Bartlett Periodogram Comparison');
legend('pwelch', 'my_{PSD}', 'Location', 'Best');

figure
plot(f2, 10*log10(abs(Sy_simple)), 'LineWidth', 1.2)
hold on
grid on
plot(f2, 10*log10(abs(Sy_mySimple)), 'LineWidth', 1.2)
xlabel('Frequency (Hz)');
ylabel('10log_{10}|S_y(f)|');
title('Simple Periodogram Comparison')
legend('pwelch', 'my_{PSD}', 'Location', 'Best');

function [Sx, f] = my_PSD(x, N, varargin)
    if nargin == 3
        M = varargin{1};
        if mod(N,M) ~= 0
            L = ceil(N/M)*M - N;
            x = [x, zeros(1, L)];
            N = N + L;
        end
        D = N/M;
        Sx = zeros(M, D);
        f = -1/2:1/D:1/2-1/D;
        for i = 1:M
            Sx(i, :) = 1/D*abs(my_DFT(x(D*(i - 1) + 1: D*i), length(x(D*(i - 1) + 1: D*i)))).^2;
        end
    else
        Sx = 1/N*abs(fftshift(fft(x))).^2;
        f = -1/2:1/N:1/2-1/N;
    end
end

function X = my_DFT(x, N)
    C = 0:(N-1);
    Uf = exp((2*pi*1j*C'.*C)/N);
    X = Uf*x(:);
    X = circshift(X, floor(N/2));
end