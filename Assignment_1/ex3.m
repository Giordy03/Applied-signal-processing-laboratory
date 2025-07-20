clc, clear, close all;

B1 = 16;
B2 = 6;
f0 = 10;
fc = 64;
tc = 1/fc;
N = 128;
n = -N/2:N/2-1;
t = n*tc;

x = 2*B1*sinc(2*B1*t) + (4*B2*sinc(2*B2*t) - 3*B2*(sinc(B2*t)).^2).*cos(2*pi*f0*t);

T0 = N*tc;
Df = fc/N;

M = 100*N;
Df_new = fc/M;
f = (-M/2:M/2-1)*Df_new;
Y_f = rectangularPulse(-B1, B1, f) + rectangularPulse(f0 - B2, f0 + B2, f) + rectangularPulse(- f0 ...
    - B2, - f0 + B2, f) - 3/2*(triangularPulse(f0-B2, f0, f0+B2, f) + triangularPulse(-f0-B2, -f0, -f0+B2, f));

%compute dft
k = (-N/2:N/2-1)*Df;
DFT = fftshift(fft(x, N));

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k, tc*abs(DFT), 'b')
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT", N))
legend("|Y(f)|", "|X(k)|")

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k, tc*abs(DFT), 'b')
xlim([-5 5])
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT, x-axis limited", N))
legend("|Y(f)|", "|X(k)|")

x_a = [x, zeros(1, 128)];
x_b = [x, zeros(1, 384)];

N_a = length(x_a);
N_b = length(x_b);

Df_a = fc/N_a;
Df_b = fc/N_b;

k_a = (-N_a/2:N_a/2-1)*Df_a;
k_b = (-N_b/2:N_b/2-1)*Df_b;

%compute dft
DFT_a = fftshift(fft(x_a, N_a));
DFT_b = fftshift(fft(x_b, N_b));

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k_a, tc*abs(DFT_a), 'b')
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT", N_a))
legend("|Y(f)|", "|X_a(k)|")

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k_a, tc*abs(DFT_a), 'b')
xlim([-5 5])
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT, x-axis limited", N_a))
legend("|Y(f)|", "|X(k)|")

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k_b, tc*abs(DFT_b), 'b')
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT", N_b))
legend("|Y(f)|", "|X_b(k)|")

figure
plot(f, abs(Y_f), '--r', 'LineWidth', 1.5)
grid on; hold on
stem(k_b, tc*abs(DFT_b), 'b')
xlim([-5 5])
xlabel("f [Hz]"), ylabel("Magnitude")
title(sprintf("Frequency domain analysis, using %d points for DFT, x-axis limited", N_b))
legend("|Y(f)|", "|X(k)|")

w_t = rectangularPulse(-T0/2, T0/2, t);
W_f = fftshift(fft(w_t));
Z_f = conv(Y_f, W_f);

% find first and last non zero element in Y_f
Y_f_nonZero = find(Y_f ~= 0);
first = Y_f_nonZero(1);
last = Y_f_nonZero(end);
Bw = f(last)-f(first);
f_tilde = f(first:last);
Y_tilde = Y_f(first:last);

Q = M*(1+Bw/fc);
Df_q = (fc+Bw)/Q;
f_q = -(fc+Bw)/2:Df_q:(fc+Bw)/2-(fc+Bw)/Q;

W_f = T0*sinc(f_q*T0);

Z_f_trapz = zeros(1, length(f_tilde));
for i = 1:length(f)
    result = Y_tilde .* W_f(i:length(Y_tilde)+i-1);
    Z_f_trapz(i) = trapz(f_tilde, result);
end

figure
hold on, grid on
plot(f, abs(Z_f_trapz), "LineWidth", 2)
stem(k, tc*abs(DFT))
xlabel("f [Hz]"), ylabel('Magnitude')
title("Convolution and DFT")
legend('|Z(f)|','|X(k)|')

figure
hold on, grid on
plot(f, abs(Z_f_trapz), "LineWidth", 2)
stem(k, tc*abs(DFT))
xlim([-5 5])
xlabel("f [Hz]"), ylabel('Magnitude')
title("Convolution and DFT, x-axis limited")
legend('|Z(f)|','|X(k)|')

