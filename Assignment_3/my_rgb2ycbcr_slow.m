function YCbCrimg = my_rgb2ycbcr_slow(RGBimg)
    YCbCrimg = zeros(size(RGBimg));
    [M, N, ~] = size(RGBimg);
    matrix = [0.299, 0.587, 0.114; -0.169, -0.331, 0.5; 0.5, -0.419, -0.081];
    for i = 1:M
        for j = 1:N
            YCbCrimg(i, j, :) = matrix*squeeze(RGBimg(i,j,:));
        end
    end
end 