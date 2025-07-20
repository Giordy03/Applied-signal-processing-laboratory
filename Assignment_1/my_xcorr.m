function z = my_xcorr(x, y)
    length_xcorr = length(x) + length(y) - 1;
    
    %because we want always row vectors
    x = x(:)';
    y = y(:)';

    x_padded = [zeros(1, length(y)-1), x, zeros(1, length(y)-1)];

    z = zeros(1, length_xcorr);

    for n = 1:length_xcorr
        x_k = x_padded(n:n + length(y) - 1);
        z(n) = conj(x_k) * y';
    end
end