
clc
close all
clear 

Tc = 2;
T_tot = 5;
Tb = 1e-3;

fc = 4e3;
tc = 1/fc;
N = T_tot*fc;
t = (0:N-1)*tc;

B = floor(Tc/Tb);
k = 1:B;
b = binornd(1, 0.5, 1, B);

%map zeros to -1
index0 = find(b==0);
b(index0) = -1;

c = repmat(b, Tb*fc, 1); %??? rivedi
c = c(:);

M = length(c);
Td = round((T_tot - Tc)*rand(1), 3);
fprintf('The generated delay is %.3f s\n', Td)
s = [zeros(Td*fc, 1); c; zeros(int32((T_tot-Tc-Td)*fc), 1)];

dev = sqrt(50);
w = dev*randn(N,1);

x = w + s;

figure
sgtitle("Time analysis")
subplot(2,1,1)
plot(t, s)
xlabel("t [s]"), ylabel("s(t)")
subplot(2,1,2)
plot(t, x)
xlabel("t [s]"), ylabel("x(t)")

%cross correlation
z_practical = my_xcorr(x, c); 

time_axis = (-M+1:N-1)*tc;

figure
plot(time_axis, z_practical, "LineWidth", 2)
grid on
xlabel("t [s]"), ylabel("Amplitude")
title("my_xcorr(x, c), noise variance=50")

%find the delay
[~, idx] = max(z_practical);
T_d = time_axis(idx); 
fprintf('Variance: 50 -> The found delay is %.3f s\n', T_d)

%% rerun the script changing the standard deviation
dev = sqrt(100);
w = dev*randn(N,1);

x = w + s;

figure
sgtitle("Time analysis, noise variance = 100")
subplot(2,1,1)
plot(t, s)
xlabel("t [s]"), ylabel("s(t)")
subplot(2,1,2)
plot(t, x)
xlabel("t [s]"), ylabel("x(t)")

%cross correlation
z_practical = my_xcorr(x, c); 

z_theory = xcorr(x, c); 
time_axis = (-M+1:N-1)*tc;

figure
plot(time_axis, z_practical, "LineWidth", 2)
xlabel("time"), ylabel("Amplitude")
title(sprintf('My_xcorr, noise variance = %d', dev^2))
grid on

% find the delay
[d_hat, idx] = max(z_practical);
T_d = time_axis(idx); 

fprintf('Variance: 100 -> The found delay is %.3f s\n', T_d)


%% PART 3
P = 20;
f0 = randi(P)*10;

var_value = 50;
dev = sqrt(var_value);
w = dev*randn(N,1);

fprintf('\nVariance = %d: \n', var_value);
fprintf('Randomly generated value of frequency f_0: %.2f Hz\n', f0);
fprintf('Randomly generated value of delay T_D: %.3f s\n', Td);

t = (0:length(c)-1)*tc; 
p = c.*cos(2*pi*f0*t');

s = [zeros(Td*fc, 1); p; zeros(int32((T_tot-Tc-Td)*fc), 1)];

y = s + w;
figure
plot(y)
f0_vector = (1:P)*10;

A = zeros(P, N + M - 1);

for i=1:P
    f0 = f0_vector(i);
    p = c.*cos(2*pi*f0*t');
    A(i, :) = my_xcorr(y, p);
end

lags = (-M+1:N-1)*tc;

figure
mesh(lags, f0_vector, A);
xlabel('Lag [s]');
ylabel('Frequency [Hz]');
zlabel('Cross-correlation');
title(sprintf('Cross-correlation Surface (Variance = %d)', var_value));
grid on

[max_col, d] = max(A, [], 2); 
[~, f0_idx] = max(max_col);

delay_max = lags(d(f0_idx)-1); % -1 because matlab has indices starting from 1??
f0_max = f0_vector(f0_idx);

fprintf('Estimated Delay of Max: %.3f s\n', delay_max);
fprintf('Estimated Frequency of Max: %.2f Hz\n', f0_max);

%% rerun the script changing the standard deviation
var_value = 100;
dev = sqrt(var_value);

P = 20;
f0 = randi(P)*10;
fprintf('\nVariance = %d: \n', var_value);
fprintf('Randomly generated value of frequency f_0: %.2f Hz\n', f0);
fprintf('Randomly generated value of delay T_D: %.3f s\n', Td);

t = (0:length(c)-1)*tc; 
p = c.*cos(2*pi*f0*t');

s = [zeros(Td*fc, 1); p; zeros(int32((T_tot-Tc-Td)*fc), 1)];

w = dev*randn(N,1);

y = s + w;

f0_vector = (1:P)*10;

A = zeros(P, N + M - 1);

for i=1:P
    f0 = f0_vector(i);
    p = c.*cos(2*pi*f0*t');
    A(i, :) = my_xcorr(y, p);
end

lags = (-M+1:N-1)*tc;

figure
mesh(lags, f0_vector, A);
xlabel('Lag [s]');
ylabel('Frequency [Hz]');
zlabel('Cross-correlation');
title(sprintf('Cross-correlation Surface (Variance = %d)', var_value));
grid on

[max_col, d] = max(A, [], 2); 
[~, f0_idx] = max(max_col);

delay_max = lags(d(f0_idx)-1); % -1 because matlab has indices starting from 1??
f0_max = f0_vector(f0_idx);

fprintf('Estimated Delay of Max: %.3f s\n', delay_max);
fprintf('Estimated Frequency of Max: %.2f Hz\n', f0_max);

