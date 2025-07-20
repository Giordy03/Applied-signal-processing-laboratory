clc; clear; close all;

x = randi(50, [1, 20]);
y = randi(30, [1, 30]);

max_xy = max(length(x), length(y));
t_circ = 0:max_xy-1;

%% Circular Convolution
z = my_conv(x, y, 'c');
z_theory = cconv(x, y, max_xy);
n = (0:length(z)-1);

figure;
hold on;
grid on;
stem(n, z_theory, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
stem(n, z, 'o--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Samples')
ylabel('Amplitude')
title('(b)')
legend('Theoretical with cconv()', 'Obtained with my\_conv()')
title('Circular Convolution')
%saveas(gcf, 'dataes1/circular_convolution.png');


%% Linear convolution
z = my_conv(x,y, 'l');
z_theory = conv(x,y);
n = (0:length(z)-1);

figure;
hold on;
grid on;
stem(n, z_theory, 'LineWidth', 1.5, 'MarkerSize', 6);
stem(n, z, 'o--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('samples')
ylabel('Amplitude')
title('(a)')
legend('Theoretical with cconv()', 'Obtained with my\_conv()')
title('Linear Convolution')
%saveas(gcf, 'dataes1/linear_convolution.png');


%check linear convolution also with a lower number of inputs in the
%function
z = my_conv(x,y);

figure;
hold on;
grid on;
stem(n, z_theory, 'LineWidth', 1.5, 'MarkerSize', 6);
stem(n, z, 'o--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('time (s)')
ylabel('Amplitude')
title('(a)')
legend('Theoretical with cconv()', 'Obtained with my\_conv()')
title('Linear Convolution (checking the function my\_conv() with one input less)')

%% Generate triangular sequence
N = 45; 
M = (N+1)/2;

% two rectangular sequences used to generate the triangle
rect1 = 1/sqrt(M)*ones(1, M);
rect2 = rect1;
x = my_conv(rect1, rect2);
n = 0:N-1;

figure;
stem(n, x, 'b', 'LineWidth', 1.5);
hold on;
grid on;
xlabel('n');
ylabel('Amplitude');
title('Triangular signal obtained with my\_conv()');
%saveas(gcf, 'dataes1/triangular_sequence.png');


%% DFT of the triangular sequence

X_k = fftshift(fft(x));
f = -1/2 + 1/(2*N):1/N:1/2 - 1/(2*N);

figure;
stem(f, abs(X_k), 'b', 'LineWidth', 1.5);
hold on;
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('|X(k)| with X(k) as the DFT of a triangular signal');
%saveas(gcf, sprintf('dataes1/DTFT_triangular.png'));


% adding zeros to the triangular sequence and comparing it with the DTFT 
[X_k1, DTFT1] = compute_DTFT(64, x, M);
[X_k2, DTFT2] = compute_DTFT(128, x, M);
[X_k3, DTFT3] = compute_DTFT(256, x, M);

function [X_k, DTFT] = compute_DTFT(N, x, M)
    x = [x; zeros(N - length(x), 1)];
    X_k = fftshift(fft(x));
    f = -1/2:1/N:1/2 - 1/(N);
    f_DTFT = linspace(-1/2, 1/2-1/N, 4000);
    DTFT = (sin(pi*f_DTFT*M).^2)./(M*(sin(pi*f_DTFT).^2));

    figure
    stem(f, abs(X_k));
    hold on
    grid on
    plot(f_DTFT, DTFT, 'LineWidth', 1.5);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    legend('DFT |X_k|', 'Theoretical DTFT') 
    title(sprintf('N = %d', N))
    %saveas(gcf, sprintf('dataes1/DTFT_N%d.png', N));
    title(sprintf('Comparison between DFT and Theoretical DTFT - N = %d', N))
    
    %show zoomed version
    figure
    stem(f, abs(X_k));
    hold on
    grid on
    plot(f_DTFT, DTFT, 'LineWidth', 1.5);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    legend('DFT |X_k|', 'Theoretical DTFT')
    xlim([-0.1, 0.1])
    ylim([0, 6])
    title(sprintf('N = %d', N))
    %saveas(gcf, sprintf('dataes1/zoom_DTFT_N%d.png', N));
    title(sprintf('Comparison between DFT and Theoretical DTFT - N = %d (zoomed)', N)) 
end
