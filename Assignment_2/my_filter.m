function y = my_filter(b, a, x)
    if any(abs((roots(a))) >= 1)
        errordlg('The filter is not stable!', 'Error');
        y = [];
        return
    end
    Na = length(a);
    Nb = length(b);
    l_x = length(x);
    max_l = max(Na, Nb);
    x = [zeros(1, max_l - 1), x(:)'];
    y = zeros(1, length(x));

   
    b_flipped = fliplr(b);
    a_flipped = fliplr(a);
    for i = 1:length(x)
        y(i+Na-1) = sum(b_flipped(1:min(end, length(x)-i+1)).*x(i:min(Nb+i-1, end))) - sum(a_flipped(1:min(end-1, length(y) - i)).*y(i: min(Na+i-2, end-1)));
    end
    y = y(max_l:max_l + l_x - 1);
end
