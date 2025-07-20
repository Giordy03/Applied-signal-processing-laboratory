clc; clear; close all;

syms t f B
assume(in(t, 'real'));
assume(in(B, 'real') & B>0); % B is positive
assume(in(f, 'real'))

x_t = 32*t*exp(-8*t)*heaviside(t);
% check the energy
En = int(abs(x_t)^2, t, 0, inf);
fprintf("The energy of the signal is %.3f.\n", En)

% obtain Fourier Transform
sympref('FourierParameters', [1, -2*pi]);
X_f = fourier(x_t, f);
disp(['The Fourier Transform is: ', char(X_f)]);


% check bandwidth | 99.9% of energy inside
EnX = int(abs(X_f)^2, f, -B, B); 
s = vpasolve(EnX == 0.999*En, B);
Bx = eval(s);
fprintf("The bandwidth from 99.9%% of energy is %.3f Hz.\n", Bx)
fs = 6*Bx;

% find T0
c = solve(x_t<1e-6, t);
T0 = eval(c);
fprintf("The positive solution for T0 is %.3f s.\n", T0(1))

% find N, round it to the nearest power of 2, find again T0 accordingly
ts = 1/fs;
N = T0(1)/ts; % because time must be greater than 0 so we take t0 = 3.26
fprintf("N is equal to %.2f.\n", N)

N = 2^nextpow2(N);
fprintf("N is rounded to %d.\n", N)
T0 = N*ts;
Df = 1/T0;
fprintf("The updated value for T0 is %.3f s.\n", T0)
fprintf("The Delta f is %.3f Hz.\n", Df)

n = (0:N-1)*ts;
x_samp = 32*n.*exp(-8.*n);

k = (-N/2:N/2-1)*Df;
DFT = fftshift(fft(x_samp, N));

M = 100*N;
t = (0:M-1)*ts/100;
Df_new = fs/M; 
f = (-M/2:M/2-1)*Df_new;
fprintf("The Delta f with M points is: %.3f Hz.\n", Df_new)

x_t_fun = matlabFunction(x_t);
X_f_fun = matlabFunction(X_f);

figure
plot(t, x_t_fun(t), '--r', 'LineWidth', 2);
grid on; zoom on; hold on;
stem(n, x_samp, 'b');
xlabel('t [s]'); ylabel('Magnitude');
legend('x(t) continuous', 'x[n] discrete')
title('Continuous and Discrete Representation of x(t)');
%saveas(gcf, 'dataes2/x_t.png');


figure
plot(f, abs(X_f_fun(f)), '--r', 'LineWidth', 2);
grid on; zoom on; hold on;
stem(k, ts*abs(DFT), 'b');
xlabel('f [Hz]'); ylabel('Magnitude');
title('(a)')
legend('|X(f)|', 'Magnitude of the DFT |X(k)|')
title('Magnitude of X(f) and DFT |X(k)|');
%saveas(gcf, 'dataes2/dft.png');

figure
plot(f, abs(X_f_fun(f)), '--r', 'LineWidth', 2);
grid on; zoom on; hold on;
stem(k, ts*abs(DFT), 'b');
xlim([8, 26])
ylim([0, 0.03])
title('(b)')
xlabel('f [Hz]'); ylabel('Magnitude');
legend('|X(f)|', 'Magnitude of the DFT |X(k)|')
title('Magnitude of X(f) and DFT |X(k)| (zoomed version)');
%saveas(gcf, 'dataes2/zoom_dft.png');

