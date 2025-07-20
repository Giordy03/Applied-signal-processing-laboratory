function est_pitch = my_pitch_estimator(x, fs, pl_on)
    [R, lags] = xcorr(x, x, 'coeff');

    % Consider only positive lags 
    pos_idx = lags > 0;
    R = R(pos_idx);
    lags = lags(pos_idx);

    [pks, locs] = findpeaks(R, 'MinPeakHeight', 0.05);

    [~, best_idx] = max(pks); % since the pitch is the maximum peak
    best_lag = lags(locs(best_idx));
    est_pitch = fs / best_lag;
    
    if pl_on
        t_lags = lags / fs;
        figure;
        plot(t_lags, R);
        hold on;
        xline(best_lag / fs, 'r--');
        text(best_lag/fs, R(lags == best_lag), sprintf('f = %.2f Hz, t = %.d s', est_pitch, best_lag/fs), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        xlabel('Lag (s)');
        ylabel('Autocorrelation');
        title('Pitch Estimation via Autocorrelation');
        ylim([min(R), 1.1]);
        grid on;
    end
end

