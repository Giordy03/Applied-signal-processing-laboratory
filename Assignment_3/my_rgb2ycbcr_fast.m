function YCbCrimg = my_rgb2ycbcr_fast(RGBimg)
R = double(RGBimg(:, :, 1));
G = double(RGBimg(:, :, 2));
B = double(RGBimg(:, :, 3));

Y  = 0.299 * R + 0.587 * G + 0.114 * B;
Cb = -0.169 * R - 0.331 * G + 0.5   * B;
Cr = 0.5   * R - 0.419 * G - 0.081 * B;

YCbCrimg(:,:,1) = Y;
YCbCrimg(:,:,2) = Cb;
YCbCrimg(:,:,3) = Cr;
end