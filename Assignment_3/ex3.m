clear
clc
close all

parrots = imread('parrots.jpg');
parrots = double(parrots)/255; 
figure
imshow(parrots, [])

tic
parrots_slow = my_rgb2ycbcr_slow(parrots);
elapsed_time_slow = toc;
fprintf('Elapsed time slow: %d.\n', elapsed_time_slow)

tic
parrots_fast = my_rgb2ycbcr_fast(parrots); 
time_fast = toc;
fprintf('Elapsed time fast: %d.\n', time_fast)

tic
parrots_buildin = rgb2ycbcr(parrots);
time_builin = toc;
fprintf('Elapsed time build-in function: %d.\n', time_builin)

figure
sgtitle("Comparison of YCbCr components obtained using different methods")
subplot(3, 3, 1)
imshow(parrots_fast(:,:,1), [])
title("Y, fast")
subplot(3, 3, 2)
imshow(parrots_fast(:,:,2), [])
title("Cb, fast")
subplot(3, 3, 3)
imshow(parrots_fast(:,:,3), [])
title("Cr, fast")
subplot(3, 3, 4)
imshow(parrots_slow(:,:,1), [])
title("Y, slow")
subplot(3, 3, 5)
imshow(parrots_slow(:,:,2), [])
title("Cb, slow")
subplot(3, 3, 6)
imshow(parrots_slow(:,:,3), [])
title("Cr, slow")
subplot(3, 3, 7)
imshow(parrots_buildin(:,:,1), [])
title("Y, build-in")
subplot(3, 3, 8)
imshow(parrots_buildin(:,:,2), [])
title("Cb, build-in")
subplot(3, 3, 9)
imshow(parrots_buildin(:,:,3), [])
title("Cr, buld-in")

%% create filter
Q = 51;
sigma = 5;
[X, Y] = meshgrid(-Q/2:Q/2);
h = 1/(2*pi*sigma^2)*exp(-(X.^2 + Y.^2)/(2*sigma^2));
H = abs(fftshift(fft2(h, 128, 128)));
fH = -128/2:128/2-1;
figure
mesh(fH, fH, H)
title("2-D Gaussian LPF spectrum")

%% Filtering
y_parrots = parrots_buildin(:,:,1);
parrots_filtered = conv2(y_parrots, h, 'same');

figure
imshow(parrots_filtered)
title("Image filtered with 2-D Gaussian LPF")

%% Gaussian in 1 D
x = -Q/2:Q/2;
h_1d = 1/sqrt(2*pi*sigma^2)*exp(-(x.^2)/(2*sigma^2));
[N, M] = size(y_parrots);

parrots_filtered_1D = zeros(size(y_parrots));
parrots_filtered_1D_partial = zeros(size(y_parrots));

for i = 1:M
    parrots_filtered_1D_partial(:,i) = conv(y_parrots(:,i), h_1d, 'same');
end
for i = 1:N
    parrots_filtered_1D(i,:) = conv(parrots_filtered_1D_partial(i,:), h_1d, 'same');
end

figure
imshow(parrots_filtered_1D)
title("Image filtered with 1-D Gaussian LPF")
%% create the spectrum
parrot_spectrum = abs(fftshift(fft2(y_parrots)));
parrot_spectrum_filtered = abs(fftshift(fft2(parrots_filtered)));

[M, N, ~] = size(parrot_spectrum);

figure
mesh(-N/2:N/2-1, -M/2:M/2-1, 10*log10(1+abs(parrot_spectrum_filtered)))
title("Centred spectrum of the filtered image (dB scale)")

figure
mesh(-N/2:N/2-1, -M/2:M/2-1, 10*log10(1+abs(parrot_spectrum)))
title("Centred spectrum of the original image (dB scale)")

imwrite(uint8(255*parrots_filtered), 'parrots_filtered.bmp')

%% new york

new_york = imread('new_york.jpg');
new_york = double(new_york)/255;
figure
imshow(new_york, [])
title("Original RGB image")

new_york = my_rgb2ycbcr_fast(new_york);
y_new_york = new_york(:,:,1);
[M,N] = size(y_new_york);

figure
imshow(y_new_york, [])
title("Y component (my\_rgb2ycbcr\_fast)")

%% filter new york

P = 2^nextpow2(2*max(M,N));
fprintf("P = %d\n", P)
f = -P/2:P/2;
[U, V] = meshgrid(-P/2:P/2);
D = sqrt(U.^2 + V.^2);
n = 6;
D0 = P/16;

H_lp = 1./(1 + (D/D0).^(2*n));
H_hp = 1 - H_lp;


figure
mesh(f, f, H_hp)
title("6^{th} order Butterworth HPF")
%% compute fft
[M_hp, N_hp] = size(H_hp);
new_york_spectrum = fftshift(fft2(y_new_york, M_hp, N_hp)); % non sappiamo per la size
new_york_filtered_fft = new_york_spectrum.*H_hp;

new_york_filtered = real(ifft2(new_york_filtered_fft));
new_york_filtered = new_york_filtered(1:M,1:N);

figure
imshow(new_york_filtered)
title("Filtered image with HPF")

%%
[M, N, ~] = size(new_york_spectrum);
figure
mesh(-N/2:N/2-1, -M/2:M/2-1, 10*log10(1+abs(new_york_filtered_fft)))
title("Centred spectrum of filtered image")
%%
figure
mesh(-N/2:N/2-1, -M/2:M/2-1, 10*log10(1+abs(new_york_spectrum)))
title("Centred spectrum of original image")
imwrite(uint8(255*abs(new_york_filtered)), 'newyork_filtered.bmp')
