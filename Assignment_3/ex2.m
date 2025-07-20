%% Exercise 2 

clear
close all
clc

A4 = 440;
fs = 44100;
T = 2.5;
t = 0:1/fs:T-1/fs;

% Chord
C5 = A4*2^(3/12);
chord1 = chord_generator(C5, 'maj');
chord1_obj = audioplayer(chord1, fs);
playblocking(chord1_obj)

% Chord progression
A4 = 440;
chord_A4 = chord_generator(A4, 'maj');
D5 = A4*2^(5/12);
chord_D5 = chord_generator(D5, 'maj');
B4 = A4*2^(2/12);
chord_B4m = chord_generator(B4, 'min');
F4sharp = A4*2^(-3/12);
chord_F4sharpm = chord_generator(F4sharp, 'min');
G4 = A4*2^(-2/12);
chord_G4 = chord_generator(G4, 'maj');

my_progression = [chord_D5, chord_A4, chord_B4m, chord_F4sharpm, chord_G4, chord_D5, chord_G4, chord_A4];
my_progression_obj = audioplayer(my_progression, fs);
playblocking(my_progression_obj)
audiowrite("my_progression.wav", my_progression, fs);

t_prog = (0:length(my_progression)-1)/fs;
figure;
plot(t_prog, my_progression);
hold on;

% Chord length in seconds:
T = 2.5;
for k = 1:7
    xline(k*T, 'k--','LineWidth',1);
end
xlabel('Time (s)');
ylabel('Amplitude');
title('Full chord progression waveform');
grid on;

% Cathedral effects
[data, fs_cathedral] = audioread("impulse_revcathedral.wav");
r = fs/fs_cathedral; 
n = 4;
cutoff = 1/r;
imp_res_cath = interp(data, r, n, cutoff);
M = 2^(nextpow2(length(imp_res_cath)+length(my_progression)));
fft_cath = fft(imp_res_cath, M);
fft_progr = fft(my_progression, M);
fft_filt = fft_cath.'.*fft_progr;
filtered_progr = real(ifft(fft_filt,M));
filtered_progr_norm = filtered_progr/max(abs(filtered_progr));

filtr = audioplayer(filtered_progr_norm, fs);
playblocking(filtr)
audiowrite("my_progression_revcathedral.wav", filtered_progr_norm, fs);

% Time-domain comparison
t_prog = (0:length(my_progression)-1)/fs;
t_rev = (0:length(filtered_progr_norm)-1)/fs;

figure;
plot(t_prog, my_progression, 'b');
hold on;
plot(t_rev, filtered_progr_norm, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Reverberated Progression');
legend('Original','With Reverb');
grid on;

% Frequency-domain comparison
f = (0:M/2-1) * (fs/M);  
mag_orig = 10*log10(abs(fft_progr));
mag_rev  = 10*log10(abs(fft_filt));

figure;
semilogx(f, mag_orig(1:M/2), 'b', f, mag_rev(1:M/2), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Original vs Reverberated Progression Spectrum');
legend('Original','With Reverb');
grid on;

t_imp = (0:length(imp_res_cath)-1)/fs;
figure;
subplot(2,1,1)
plot(t_imp, imp_res_cath);
xlabel('Time (s)');
ylabel('Amplitude');
title('Interpolated cathedral impulse response');
grid on;
Nfft = 2^nextpow2(length(imp_res_cath));
H = fft(imp_res_cath, Nfft);
f = (0:Nfft/2-1) * (fs/Nfft);
subplot(2,1,2)
semilogx(f,10*log10(abs(H(1:Nfft/2))));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude response of cathedral impulse');
grid on;

function note = note_generation(f0)
    fs = 44100;
    T = 2.5;
    t = 0:1/fs:T-1/fs;
    lfo = f0/200;

    saw = sawtooth(2*pi*f0*t);
    sine = sin(2*pi*lfo*t);
    pulse = saw < sine;
   
    wp = f0;
    ws = 16*f0;
    Bt = 15*f0;
    fc_norm = (wp+ws)/(2*fs); 
    
    N = 6.1*pi/(Bt/fs*2*pi);
    M = ceil((N-1)/2);
    n = -M:M;

    h_id = 2*fc_norm*sinc(2*fc_norm*(n-M));
    w_bartlett = 1 - 2*abs(n)/M;
    h = h_id.*w_bartlett;
    h_norm = h/sum(h);
    note = filter(h_norm, 1, pulse);
end

function smooth_chord = chord_generator(root_note, type)
    T = 2.5;
    fs = 44100;
    t = 0:1/fs:T-1/fs;

    switch type
      case 'maj'
        f2 = 2^(4/12) * root_note;   % major third
      case 'min'
        f2 = 2^(3/12) * root_note;   % minor third
      otherwise
        error('chord_generator: invalid type, use ''maj'' or ''min''');
    end
    f3 = 2^(7/12)*root_note;

    note1 = note_generation(root_note);
    note2 = note_generation(f2);
    note3 = note_generation(f3);
    
    chord = note1 + note2 + note3;
    chord_norm = chord./max(abs(chord));

    % ADSR envelope
    ADSR = interp1([0, 0.16, 0.32, 0.6, 1]*T, [0, 1, 0.7, 0.7, 0], t, 'pchip');
    smooth_chord = chord_norm.*ADSR;

    figure;
    plot(t, chord, 'b');
    hold on;
    plot(t, ADSR * max(abs(chord)), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Chord waveform with ADSR envelope');
    legend('Chord','ADSR envelope');
    grid on;

    L    = length(chord);
    NFFT = 2^nextpow2(L);
    X    = fft(chord, NFFT);
    f    = fs/2 * linspace(0,1,NFFT/2+1);

    figure;
    semilogx(f, 20*log10(abs(X(1:NFFT/2+1))));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Magnitude spectrum of single chord');
    grid on;
end

