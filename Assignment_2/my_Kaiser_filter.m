function [Bk, N, beta]  = my_Kaiser_filter(As, Bt, fs, fc, type)

    fc_norm = fc/fs;
    N = ceil((As-8)/(2.285*(Bt/fs)*2*pi));
    M = ceil((N-1)/2); 
    n = 0:(N-1);

    if As > 50
        beta = 0.1102*(As-8.7);
    elseif As >= 21
        beta = 0.5842*(As-21)^0.4 + 0.07886*(As-21);
    else
        beta = 0;
    end

    if strcmp(type, '-lp')
         h_id = 2*fc_norm*sinc(2*fc_norm*(n-M));

    elseif strcmp(type, '-hp')
        h_id = (n==M) - 2*fc_norm*sinc(2*fc_norm*(n-M)); %n == M corresponds to delta

    elseif strcmp(type, '-bp')
         h_id = 2*fc_norm(2)*sinc(2*fc_norm(2)*(n-M)) - 2*fc_norm(1)*sinc(2*fc_norm(1)*(n-M));

    elseif strcmp(type, '-bs')
         h_id = (n==M) - (2*fc_norm(2)*sinc(2*fc_norm(2)*(n-M)) - 2*fc_norm(1)*sinc(2*fc_norm(1)*(n-M)));

    else
        errordlg('Invalid filter type', 'ERROR')
        return
    end

    %w_Kaiser = kaiser(N, beta)'; 
    w_Kaiser = besseli(0, beta * sqrt(1 - ((2*n - N + 1)/(N - 1)).^2)) / besseli(0, beta);
    Bk = h_id.*w_Kaiser;

end