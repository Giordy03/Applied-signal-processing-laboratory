clear
clc
close all

[audio_guitar, fs] = audioread("A2_guitar.wav");
audio_guitar = audio_guitar(:,1);

audio_guitar = audio_guitar/max(audio_guitar);
player = audioplayer(audio_guitar, fs);
playblocking(player);

t = 0:1/fs:length(audio_guitar)/fs -1/fs;

figure
plot(t, audio_guitar)
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Time-Domain Representation of the Guitar Signal')

% periodogram
NFFT = 2000;
N = length(audio_guitar) / NFFT;
window = hamming(NFFT);
n_overlap = NFFT/2;
[S_guitar, f_guitar] = pwelch(audio_guitar, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_guitar, 10*log10(S_guitar))
grid on
xlabel('Frequency [Hz]')
ylabel('PSD')
title('Power Spectral Density of Guitar Signal (Welch Method)')

% spectrogram
N_window = 4000;
N = length(audio_guitar) / N_window;
window = hamming(N_window);
n_overlap = N_window*9/10;
NFFT = 2^12;

figure
spectrogram(audio_guitar, window, n_overlap, NFFT, fs, 'yaxis')
ylim([0 8])
title('Spectrogram of the Guitar Signal')


est_pitch_guitar = my_pitch_estimator(audio_guitar, fs, 1);




%% file audio 2
[audio_bending, fs] = audioread("Bending.wav");
audio_bending = audio_bending(:,1);

audio_bending = audio_bending/max(audio_bending);
player = audioplayer(audio_bending, fs);
playblocking(player);

t = 0:1/fs:length(audio_bending)/fs -1/fs;

figure
plot(t, audio_bending)
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Time-Domain Representation of the Bending signal')


% periodogram
NFFT = 2000;
N = length(audio_bending) / NFFT;
window = hamming(NFFT);
n_overlap = NFFT/2;
[S_guitar, f_guitar] = pwelch(audio_bending, window, n_overlap, NFFT, fs, 'centered');

figure
plot(f_guitar, 10*log10(S_guitar))
grid on
xlabel('Frequency [Hz]')
ylabel('PSD')
title('Power Spectral Density of Bending Signal (Welch Method)')

% spectrogram
N_window = 4000;
N = length(audio_bending) / N_window;
window = hamming(N_window);
n_overlap = N_window*9/10;
NFFT = 2^12;
note1 = audio_bending(0.22*fs:1.26*fs, 1);
note1_pitch = my_pitch_estimator(note1, fs, true);
note2 = audio_bending(1.26*fs:2.1*fs, 1);
note2_pitch = my_pitch_estimator(note2, fs, true);
note3 = audio_bending(2.1*fs:3.6*fs, 1);
figure
spectrogram(audio_bending, window, n_overlap, NFFT, fs, 'yaxis')
ylim([0 8])
title('Spectrogram of the Bending Signal')



N = length(note3);               
M = round(0.1 * fs);   
L = round(0.9 * M);                       
S = floor((N - L)/(M - L));       

pitches_ind = zeros(2,S);
for i = 1:S
    start_idx = (i-1)*(M-L) + 1;
    end_idx   = start_idx + M - 1;
    frame = note3(start_idx:end_idx);
    pitches_ind(1,i) = my_pitch_estimator(frame, fs, false);
    pitches_ind(2,i) = start_idx + floor((length(frame)-1)/2); 
end

figure
plot(pitches_ind(2,:), pitches_ind(1,:)) 
grid on
xlabel('Time (s)')
ylabel('Estimated pitch (Hz)')
title('Pitch estimation over time for bent note')

fc = 5e3;
% window = hamming

%from table -> min attenuation: 50 dB -> Hamming
trans_band = 500;
trans_band_norm = trans_band/(fs);

N = ceil(6.6*pi/(2*pi*trans_band_norm));
fprintf('The number of points is %d.\n', N)

fc_norm = fc/(fs); 

fprintf('The normalized (with respect to f_s) cutoff frequency is %.2f.\n\n', fc_norm)

n = 0:(N-1);
w_Hamming = 0.54 - 0.46*cos(2*pi*n/(N-1));
df = fs/N;
f = -fs/2:df:fs/2-df;

M = ceil((N-1)/2);
h_id = 2*fc_norm*sinc(2*fc_norm*(n-M));

h = h_id.*w_Hamming;

imp_res = impz(h,1,1024, fs);
[H, ~] = freqz(imp_res, 1, 1024, fs);
n_imp = 0:length(imp_res)-1;
bending_filtered = filter(h, 1, audio_bending);

% undersample -> factor 4
under_factor = 4;
fs_2 = fs/under_factor;

bending_filtered_under = bending_filtered(1:under_factor:end);

% spectrogram
N_window = 4000;
N = length(bending_filtered_under) / N_window;
window = hamming(N_window);
n_overlap = N_window*9/10;
NFFT = 2^12;

figure
spectrogram(bending_filtered_under, window, n_overlap, NFFT, fs_2, 'yaxis')
title('Spectrogram bending filter - undersample factor: 4')

% play it
audio_normalized = bending_filtered_under / max(abs(bending_filtered_under));
player = audioplayer(audio_normalized, fs_2);
playblocking(player);


% repeat step 12
fc = 2.5e3;
% window = hamming
trans_band = 500;
trans_band_norm = trans_band/(fs);
N = ceil(6.6*pi/(2*pi*trans_band_norm));
%fprintf('Undersample factor: 4')
%fprintf('The number of points is %d.\n', N)

fc_norm = fc/(fs);

%fprintf('The normalized (with respect to f_s) cutoff frequency is %.2f.\n\n', fc_norm)

n = 0:(N-1);
w_Hamming = 0.54 - 0.46*cos(2*pi*n/(N-1));
df = fs/N;
f = -fs/2:df:fs/2-df;

M = ceil((N-1)/2);
h_id = 2*fc_norm*sinc(2*fc_norm*(n-M));
h = h_id.*w_Hamming;

imp_res = impz(h,1,1024, fs);
[H, ~] = freqz(imp_res, 1, 1024, fs);
n_imp = 0:length(imp_res)-1;
bending_filtered = filter(h, 1, audio_bending);

% undersample -> factor 8
under_factor = 8;
fs_2 = fs/under_factor;

bending_filtered_under = bending_filtered(1:under_factor:end);

% spectrogram
N_window = 4000;
N = length(bending_filtered_under) / N_window;
window = hamming(N_window);
n_overlap = N_window*9/10;
NFFT = 2^12;

figure
spectrogram(bending_filtered_under, window, n_overlap, NFFT, fs_2, 'yaxis')
title('Spectrogram bending filter - undersample factor: 8')

% play it
audio_normalized = bending_filtered_under / max(abs(bending_filtered_under));
player = audioplayer(audio_normalized, fs_2);
playblocking(player);

%% register the voice and estimate pitch

fs_voice = 8e3;
n_bits_voice = 16;

if isfile("my_speech.wav")
    my_Speech = audioread("my_speech.wav");
else
    r = audiorecorder(fs_voice, n_bits_voice, 1);
    disp('Start speaking now!')
    recordblocking(r, 5); %record for 5 seconds
    disp('End of recording.')
    y = getaudiodata(r);
    audiowrite("my_speech.wav", y, fs_voice);
    my_Speech = audioread("my_speech.wav");
end

my_Speech = my_Speech/max(abs(my_Speech));
player = audioplayer(my_Speech, fs_voice);
playblocking(player)

t = 0:1/fs_voice:length(my_Speech)/fs_voice -1/fs_voice;

figure
plot(t, my_Speech)
grid on
xlabel('Time [s]')
ylabel('Amplitude')
title('Time-Domain Representation of the Voice Signal')

% periodogram
NFFT = 2000;
N = length(my_Speech) / NFFT;
window = hamming(NFFT);
n_overlap = NFFT/2;
[S_speech, f_speech] = pwelch(my_Speech, window, n_overlap, NFFT, fs_voice, 'centered');

figure
plot(f_speech, 10*log10(S_speech))
grid on
xlabel('Frequency [Hz]')
ylabel('PSD')
title('Power Spectral Density of Voice Signal (Welch Method)')

% find only vocalized signal
my_Speech_vocalized = my_Speech(round(0.5*fs_voice)+1:round(fs_voice*1));

%estimated pitch
est_pitch_voice = my_pitch_estimator(my_Speech_vocalized, fs_voice, 1);
fprintf('The estimated pitch for the voice is: %.2f.\n', est_pitch_voice) 
