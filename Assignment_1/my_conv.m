function z = my_conv(x, y, varargin)
if (nargin == 3)
    if (varargin{1} == 'c')
    max_dim = max(length(x), length(y));
    if length(x) > length(y)
        padded = [y, zeros(1, max_dim - length(y))]';
        largest = x;
    else
        padded = [x, zeros(1, max_dim - length(x))]';
        largest = y;
    end
    M_circ = zeros(max_dim, max_dim);
    n = 0:max_dim-1;
    new_col = @(i) circshift(padded, n(i));
    for i = 1:max_dim
        M_circ(:,i) = new_col(i);
    end
    z = M_circ * largest';
    else
        length_conv = length(x) + length(y) - 1;
        x_padded = [x, zeros(1, length_conv - length(x))];
        y_padded = [y, zeros(1, length_conv - length(y))];
        z = zeros(length_conv, 1);
        for n = 1:length_conv
            x_k = x_padded(1:n);
            y_k_n = y_padded(n:-1:1);
            z(n) = x_k * y_k_n';
        end
    end
else
    length_conv = length(x) + length(y) - 1;
        x_padded = [x, zeros(1, length_conv - 1)];
        y_padded = [y, zeros(1, length_conv - 1)];
        z = zeros(length_conv, 1);
        for n = 1:length_conv
            x_k = x_padded(1:n);
            y_k_n = y_padded(n:-1:1);
            z(n) = x_k * y_k_n';
        end
end